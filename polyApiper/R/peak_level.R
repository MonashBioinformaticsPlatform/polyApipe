#
# Using peak-level data
#


se_reassign <- function(se, organism) {
    se <- load_banquet(se)
    organism <- load_banquet(organism)

    rd <- rowData(se)
    ranges <- GRanges(
            seqnames=rd$pchr,
            ranges=IRanges(start=rd$pstart, end=rd$pend),
            strand=rd$pstrand) %>%
        anchor_3p() %>%
        mutate(width=1) %>%
        join_overlap_left_directed(organism$regions)

    assert_that(length(ranges) == nrow(se))

    rowData(se)$gene_id <- ranges$gene_id
    rowData(se)$symbol <- ranges$symbol
    rowData(se)$biotype <- ranges$biotype
    rowData(se)$region <- ranges$region
    se
}
