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

    ah[[1]]
}

get_regions <- function(db, hard_extension=20, extension=2000) {
    hard_extend <- function(gr)
        gr %>% anchor_5p() %>% mutate(width=width+hard_extension) %>% trim()

    just_mcols <- function(gr)
        gr %>% mcols() %>% as.data.frame()

    # Flatten a GRangesList
    delist <- function(grl, id_col) {
        gr <- unlist(grl)
        mcols(gr)[[id_col]] <- names(gr)
        names(gr) <- NULL
        gr
    }

    left_join_mcols <- function(left, right, by) {
        df <- just_mcols(left) %>%
            mutate(.row. = row_number()) %>%
            left_join(just_mcols(right), by=by)

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
        plyranges::select(gene_id, symbol, biotype=gene_biotype, gene_strand)

    db_trans <- transcriptsBy(db) %>%
        delist("gene_id") %>%
        mutate(support=ifelse(is.na(tx_support_level),6,tx_support_level) +
                   10*(tx_biotype!="protein_coding")) %>%
        select(gene_id, tx_id, support) %>%
        left_join_mcols(db_genes, "gene_id") %>%
        filter(strand == gene_strand) %>% #Forbid antisense isoforms
        select(tx_id, gene_id, symbol, biotype, support)

    db_utr3s <- threeUTRsByTranscript(db) %>%
        delist("tx_id") %>%
        plyranges::select(tx_id) %>%
        left_join_mcols(db_trans, "tx_id")

    db_exons <- exonsBy(db) %>% delist("tx_id") %>%
        plyranges::select(tx_id) %>%
        left_join_mcols(db_trans, "tx_id")

    db_gaps <- setdiff_ranges_directed(seq_stranded_ranges, db_trans)

    db_extensions <- db_gaps %>%
        join_overlap_inner_directed(flank_downstream(db_trans,1)) %>%
        anchor_5p() %>%
        mutate(width=pmin(width, extension))



    everything <- bind_ranges(
        mutate(db_utr3s, region="3'UTR") %>% hard_extend,
        mutate(db_exons, region="exon", support=support+100) %>% hard_extend,
        mutate(db_trans, region="intron", support=support+200),
        mutate(db_extensions, region="extension", support=support+300))

    disjoinment <- disjoin_ranges_directed(everything)
    disjoinment$id <- seq_len(length(disjoinment))

    assignment <- everything %>%
        join_overlap_inner_directed(disjoinment) %>%
        just_mcols() %>%
        group_by(id) %>%
        filter(support == min(support)) %>%
        ungroup() %>%
        distinct(id, region, symbol, gene_id, biotype) %>%
        group_by(id) %>%
        mutate(good=n()==1) %>%
        ungroup()

    good <- filter(assignment, good)
    good_ranges <- disjoinment[good$id,]
    mcols(good_ranges) <- select(good, -id, -good)

    colors <- c("3'UTR"="#008800", "exon"="#000088", "intron"="#888888", "extension"="#880000")
    mcols(good_ranges)$color <- colors[good_ranges$region]

    #Fails:
    #good_ranges <- good_ranges %>%
    #    group_by(region, symbol, gene_id, biotype) %>%
    #    reduce_ranges_directed()

    write_gff3(good_ranges, "output/regions.gff", index=TRUE)


    bad <- dplyr::filter(assignment, !good)
    bad_ranges <- disjoinment[unique(bad$id),] %>% disjoin_ranges_directed()
    mcols(bad_ranges)$id <- seq_len(length(bad_ranges))
    write_gff3(bad_ranges, "output/bad.gff", index=TRUE)
}
