#
# Automated analysis steps, producing banquetry
#
# do_ functions create a file or directory rather than returning a value.
# Their inputs can also be the path to a file or directory.
#


#' @export
se_counts_weitrix <- function(se, min_reads=50, min_prop=0.05, do_vooma=FALSE, plot_filename=NULL) {
    se <- load_banquet(se)
    message(nrow(se), " x ", ncol(se), " counts")

    min_reads <- ceiling(max(min_reads, ncol(se)*min_prop))
    message("Minimum reads per gene allowed ",min_reads)

    sums <- rowSums(counts(se))
    keep <- sums >= min_reads
    se_sub <- se[keep,]

    message(nrow(se), " genes kept")

    se_sub <- computeSumFactors(se_sub)
    se_sub <- normalize(se_sub)

    if (do_vooma) {
        if (!is.null(plot_filename)) png(plot_filename)
        # TODO: better design matrix
        # TODO: low-memory solution
        result <- as_weitrix(vooma(assay(se_sub,"logcounts"), plot=TRUE, span=0.25))
        rowData(result) <- rowData(se_sub)
        colData(result) <- colData(se_sub)
        if (!is.null(plot_filename)) dev.off()
        return(result)
    } else {
        assay(se_sub, "weights") <- matrix(1,nrow=nrow(se_sub),ncol=ncol(se_sub))
        return(bless_weitrix(se_sub, "logcounts", "weights"))
    }
}

#' @export
do_se_counts_weitrix <- function(out_path, ...) working_in(out_path, {
    # High memory
    result <- se_counts_weitrix(plot_filename=file.path(out_path,"vooma.png"), ...)

    saveHDF5SummarizedExperiment(result, file.path(out_path,"weitrix"), replace=TRUE)
    invisible()
})


amalgamate <- function(vec) paste(unique(vec),collapse="/")

#' Shift weitrix from peak counts
#'
#' Convert a SummarizedExperiment of peak counts into a weitrix of shifts.
#'
#' Peaks with less than min_present non-zero counts are discarded.
#'
#' Only genes with two or more peaks are included.
#'
#' Cells with less than min_present non-missing shifts are discarded.
#'
#' @export
peaks_shift_weitrix <- function(
        peaks_se, organism,
        n_per_row=100, n_per_col=100,
        remove_mispriming=TRUE, utr_or_extension_only=TRUE) {
    peaks_se <- load_banquet(peaks_se)
    organism <- load_banquet(organism)
    message(nrow(peaks_se), " x ", ncol(peaks_se), " counts")

    peaks_se <- se_reassign(peaks_se, organism)
    n_present <- rowSums(assay(peaks_se,"counts") > 0)
    # This seems to be necessary:
    #keep <- rep(TRUE, nrow(peaks_se))
    keep <- n_present >= n_per_row
    if (remove_mispriming)
        keep <- keep & (rowData(peaks_se)$misprime == "False")
    if (utr_or_extension_only)
        keep <- keep & (
            rowData(peaks_se)$region == "3'UTR" |
            rowData(peaks_se)$region == "extension")

    keep[is.na(keep)] <- FALSE

    rd <- rowData(peaks_se)[keep,,drop=F] %>%
        as.data.frame() %>%
        rownames_to_column("name") %>%
        arrange(ifelse(pstrand == "-", -pstart, pend)) %>%
        group_by(gene_id) %>%
        filter(length(name) >= 2) %>%
        ungroup()

    grouping <- rd %>%
        select(group=gene_id, name=name, region=region)

    gene_info <- rd %>%
        group_by(gene_id) %>%
        summarize(
            symbol=amalgamate(symbol),
            biotype=amalgamate(biotype),
            regions=paste(region,collapse="-"))

    message("Kept ", nrow(rd), " peaks")

    result <- counts_shift(assay(peaks_se,"counts"), grouping)

    result <- weitrix_filter(result, n_per_row, n_per_col)

    gene_info <-
        gene_info[match(rownames(result), gene_info$gene_id),,drop=F]
    rowData(result)$symbol <- gene_info$symbol
    rowData(result)$regions <- gene_info$regions
    rowData(result)$biotype <- gene_info$biotype

    colData(result) <-
        colData(peaks_se)[match(colnames(result), colnames(peaks_se)),,drop=F]

    result
}

#' @export
do_peaks_shift_weitrix <- function(out_path, ...) working_in(out_path, {
    result <- peaks_shift_weitrix(...)
    saveHDF5SummarizedExperiment(result, file.path(out_path,"weitrix"), replace=TRUE)
    invisible()
})


#' @export
do_weitrix_components <- function(out_path, weitrix, p=20, ...) working_in(out_path, {
    weitrix <- load_banquet(weitrix)

    comp_seq <- weitrix_components_seq(weitrix, p=p, ...)
    saveRDS(comp_seq, file.path(out_path,"comp_seq.rds"))

    rand_weitrix <- weitrix_randomize(weitrix)
    rand_comp_seq <- weitrix_components_seq(rand_weitrix, p=min(p,5), ...)
    saveRDS(rand_comp_seq, file.path(out_path,"rand_comp_seq.rds"))

    invisible()
})






