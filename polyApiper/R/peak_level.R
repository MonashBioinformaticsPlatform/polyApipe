#
# Using peak-level data
#


# Union-Find algorithms to identify connected groups
unionize <- function(as,bs, further_nodes=NULL) {
    all_names <- unique(c(as,bs,further_nodes))
    as <- match(as, all_names)
    bs <- match(bs, all_names)

    group <- seq_along(all_names)
    top <- function(i) {
        if (group[i] != i)
            group[i] <<- top(group[i])
        group[i]
    }

    for(i in seq_along(as)) {
        group[top(as[i])] <- top(bs[i])
    }

    for(i in seq_along(all_names)) top(i)

    group <- match(group, unique(group))
    names(group) <- all_names
    group
}

# Helper function for se_reassign.
# If gene is protein coding, group regions based on shared transcript ids.
utr_grouping <- function(biotype, region, tx_id, misprime, orderer) {
    grouping <- rep(NA_integer_, length(region))
    if (any(is.na(biotype) | biotype != "protein_coding"))
        return(grouping)

    use <- which(!misprime & region %in% c("3'UTR","extension"))
    use <- use[ order(orderer[use]) ]
    split_tx_id <- strsplit(tx_id[use]," ")
    first_tx_id <- map_chr(split_tx_id,1)

    a <- rep(first_tx_id, times=lengths(split_tx_id))
    b <- unlist(split_tx_id)
    group <- unionize(a,b)
    grouping[use] <- group[first_tx_id]

    grouping
}


#' @export
se_reassign <- function(se, organism) {
    se <- load_banquet(se)
    organism <- load_banquet(organism)

    ranges <- rowRanges(se) %>%
        select(peak, peakdepth, misprime) %>%
        anchor_3p() %>%
        mutate(width=1) %>%
        join_overlap_left_directed(organism$regions)

    # Group 3'UTRs and extensions
    ranges <- ranges %>%
        group_by(gene_id) %>%
        mutate(
            utr_group=utr_grouping(biotype, region, tx_id, misprime=misprime, orderer=ifelse(strand=="-",-start,end)),
            region_grouped=paste0(region,ifelse(is.na(utr_group),"",utr_group))) %>%
        ungroup()

    assert_that(length(ranges) == nrow(se))

    rowData(se)$gene_id <- ranges$gene_id
    rowData(se)$symbol <- ranges$symbol
    rowData(se)$symbol_unique <- ranges$symbol_unique
    rowData(se)$biotype <- ranges$biotype
    rowData(se)$region <- ranges$region
    rowData(se)$utr_group <- ranges$utr_group
    rowData(se)$region_grouped <- ranges$region_grouped

    se
}

#' Update rowData of an SCE stored on disk
#'
#' @export
do_se_reassign <- function(se_path, organism) {
    se <- load_banquet(se_path)
    se <- se_reassign(se, organism)
    quickResaveHDF5SummarizedExperiment(se)
}
