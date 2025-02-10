#
# Regions are ranked by
# - 3'UTR (or up to 20 bases downstrand of)
# - exon  (or up to 20 bases downstrand of)
# - intron
# - extension up to 4kb
# then by
# - protein coding transcript
# - other transcripts
# then by
# - transcript support level (NA -> max+1)
#
# If this is insufficent to uniquely prioritize a region to a gene, it is dropped.
#

get_ah <- function(ahid) {
    AnnotationHub()[[ahid]]
}

#' Get AnnotationHub ID for an appropriate EnsDb
#'
#' @export
get_ensdb <- function(species=NULL, version=NULL) {
    ah <- AnnotationHub()
    ah <- ah[ ah$rdataclass == "EnsDb" ]
    ah <- ah[ ah$species == species ]
    if (!is.null(version))
        ah <- ah[map_lgl(mcols(ah)$tags, ~version %in% .)]

    if (length(ah) == 0)
        stop("Not found")

    if (length(ah) > 1) {
        print(ah)
        stop("Multiple records available")
    }

    names(ah)
}

##' Get AnnotationHub ID for an ENSEMBL genome sequence
##'
##' @export
#get_dna <- function(species, version) {
#    ah <- AnnotationHub()
#    ah <- ah[ ah$rdataclass == "TwoBitFile" ]
#    ah <- ah[ ah$species == species ]
#    pattern <- paste0(
#        "^ftp://ftp\\.ensembl\\.org/pub/release-",version,
#        "/.*\\.dna_sm\\..*\\.fa\\.gz$")
#    ah <- ah[ grepl(pattern, ah$sourceurl) ]
#
#    if (length(ah) == 0)
#        stop("Not found")
#
#    if (length(ah) > 1) {
#        print(ah)
#        stop("Multiple records available")
#    }
#
#    names(ah)
#}

get_orgdb <- function(species) {
    ah <- AnnotationHub()
    ah <- ah[ ah$rdataclass == "OrgDb" ]
    ah <- ah[ ah$species == species ]

    if (length(ah) == 0)
        stop("Not found")

    if (length(ah) > 1) {
        print(ah)
        stop("Multiple records available")
    }

    names(ah)
}



get_regions <- function(db, hard_extension=20, extension=2000, unstranded_gaps=TRUE, omit_biotypes=c()) {
    hard_extend <- function(gr) {
        # mutate produces need-to-trim warning
        suppressWarnings(
            gr %>% 
            anchor_5p() %>% 
            mutate(width=width+hard_extension) %>%
            trim()
        )
    }

    just_mcols <- function(gr)
        gr %>% mcols() %>% as.data.frame()

    # Flatten a GRangesList
    delist <- function(grl, id_col) {
        gr <- unlist(grl)
        mcols(gr)[[id_col]] <- names(gr)
        names(gr) <- NULL
        gr
    }

    inner_join_mcols <- function(left, right, by) {
        df <- just_mcols(left) %>%
            mutate(.row. = row_number()) %>%
            inner_join(just_mcols(right), by=by)

        result <- left[df$.row.,]
        mcols(result) <- dplyr::select(df, -.row.)
        result
    }

    seq_ranges <- seqinfo(db) %>% as("GRanges")
    seq_stranded_ranges <- c(
        mutate(seq_ranges, strand = "+"),
        mutate(seq_ranges, strand = "-"))

    db_genes <- genes(db) %>%
        mutate(gene_strand = strand) %>%
        plyranges::select(gene_id, symbol, biotype=gene_biotype, gene_strand) %>%
        filter(!biotype %in% omit_biotypes) %>%
        mutate(symbol_unique = uniquifyFeatureNames(gene_id, symbol))

    db_trans <- transcriptsBy(db) %>%
        delist("gene_id") %>%
        mutate(support=ifelse(is.na(tx_support_level),6,tx_support_level) +
                   10*(tx_biotype!="protein_coding")) %>%
        select(gene_id, tx_id, support) %>%
        inner_join_mcols(db_genes, "gene_id") %>%
        filter(strand == gene_strand) %>% #Forbid antisense isoforms
        select(tx_id, gene_id, symbol, symbol_unique, biotype, support)

    db_cds <- cdsBy(db) %>%
        delist("tx_id") %>%
        plyranges::select(tx_id) %>%
        inner_join_mcols(db_trans, "tx_id")

    db_utr5s <- fiveUTRsByTranscript(db) %>%
        delist("tx_id") %>%
        plyranges::select(tx_id) %>%
        inner_join_mcols(db_trans, "tx_id")

    db_utr3s <- threeUTRsByTranscript(db) %>%
        delist("tx_id") %>%
        plyranges::select(tx_id) %>%
        inner_join_mcols(db_trans, "tx_id")

    db_exons <- exonsBy(db) %>% delist("tx_id") %>%
        plyranges::select(tx_id) %>%
        inner_join_mcols(db_trans, "tx_id")

    db_gaps <- setdiff_ranges_directed(seq_stranded_ranges, db_trans)
    if (unstranded_gaps) {
        db_trans_flip <- db_trans
        strand(db_trans_flip) <- ifelse(strand(db_trans_flip)=="+","-","+")
        db_gaps <- setdiff_ranges_directed(db_gaps, db_trans_flip)
    }

    db_extensions <- db_gaps %>%
        join_overlap_inner_directed(flank_downstream(db_trans,1)) %>%
        anchor_5p() %>%
        mutate(width=pmin(width, extension))


    everything <- bind_ranges(
        mutate(db_utr3s, region="3'UTR") %>% hard_extend,
        mutate(db_cds, region="CDS", support=support+100) %>% hard_extend,
        mutate(db_utr5s, region="5'UTR", support=support+200) %>% hard_extend,
        mutate(db_exons, region="exon", support=support+300) %>% hard_extend,
        mutate(db_trans, region="intron", support=support+400),
        mutate(db_extensions, region="extension", support=support+500))

    disjoinment <- disjoin_ranges_directed(everything)
    disjoinment$id <- seq_len(length(disjoinment))

    assignment <- everything %>%
        join_overlap_inner_directed(disjoinment) %>%
        just_mcols() %>%
        group_by(id) %>%
        filter(support == min(support)) %>%
        ungroup() %>%
        group_by(id, region, symbol, symbol_unique, gene_id, biotype) %>%
        summarize(
            tx_id=paste(tx_id,collapse=" "),
            .groups="drop") %>%
        group_by(id) %>%
        mutate(good=length(id)==1) %>% # n() not working for some reason
        ungroup()

    good <- filter(assignment, good)
    good_ranges <- disjoinment[good$id,]
    mcols(good_ranges) <- select(good, -id, -good)

    colors <- c("3'UTR"="#008800", "5'UTR"="#888800", "CDS"="#008888", "exon"="#000088", "intron"="#888888", "extension"="#880000")
    mcols(good_ranges)$color <- colors[good_ranges$region]

    #Fails:
    #good_ranges <- good_ranges %>%
    #    group_by(region, symbol, gene_id, biotype) %>%
    #    reduce_ranges_directed()

    #write_gff3(good_ranges, "output/regions.gff", index=TRUE)


    bad <- dplyr::filter(assignment, !good)
    bad_ranges <- disjoinment[unique(bad$id),] %>% disjoin_ranges_directed()
    mcols(bad_ranges)$id <- seq_len(length(bad_ranges))
    #write_gff3(bad_ranges, "output/bad.gff", index=TRUE)

    list(regions=good_ranges, bad=bad_ranges)
}


#' Create a polyApiper organism directory based on ENSEMBL annotation in AnnotationHub
#'
#' @export
do_ensembl_organism <- function(
       out_path, species, version, 
       hard_extension=20, extension=2000, unstranded_gaps=TRUE, omit_biotypes=c()
       ) {
    db_ahid <- get_ensdb(species, version)
    orgdb_ahid <- get_orgdb(species)

    ensure_dir(out_path)
    config <- list(
        AnnotationHub=list(
            txdb=db_ahid, 
            orgdb=orgdb_ahid))
    write_json(config, file.path(out_path, "config.json"), auto_unbox=TRUE)

    result <- get_regions(get_ah(db_ahid),
        hard_extension=hard_extension, extension=extension, 
        unstranded_gaps=unstranded_gaps, omit_biotypes=omit_biotypes)
    write_gff3(result$regions, file.path(out_path,"regions.gff3"), index=TRUE)
    write_gff3(result$bad, file.path(out_path,"bad.gff3"), index=TRUE)

    invisible()
}



##' @export
#do_organism <- function(out_path, txdb,
#        hard_extension=20, extension=2000, unstranded_gaps=TRUE, omit_biotypes=c()) {
#    ensure_dir(out_path)
#    
#    saveRDS(txdb, file.path(out_path,"txdb.rds"))
#    
#    result <- get_regions(txdb,
#        hard_extension=hard_extension, extension=extension, 
#        unstranded_gaps=unstranded_gaps, omit_biotypes=omit_biotypes)
#    
#    write_gff3(result$regions, file.path(out_path,"regions.gff3"), index=TRUE)
#    write_gff3(result$bad, file.path(out_path,"bad.gff3"), index=TRUE)
#    
#    invisible()
#}

