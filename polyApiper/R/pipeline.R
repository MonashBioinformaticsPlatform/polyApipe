#
# Automated analysis steps, producing banquetry
#
# do_ functions create a file or directory rather than returning a value.
# Their inputs can also be the path to a file or directory.
#

do_ensembl_organism <- function(
        out_path, species, version, extension=2000) {

}


#' @export
se_counts_weitrix <- function(se, min_reads=max(10, 0.1*ncol(se))) {
    se <- load_banquet(se)

    sums <- rowSums(counts(se))
    keep <- sums >= min_reads
    se_sub <- se[keep,]

    se_sub <- computeSumFactors(se_sub)
    se_sub <- normalize(se_sub)

    #TODO: weights
    assay(se_sub, "weights") <- matrix(1,nrow=nrow(se_sub),ncol=ncol(se_sub))

    bless_weitrix(se_sub, "logcounts", "weights")
}

#' @export
do_se_counts_weitrix <- function(out_path, ...) working_in(out_path, {
    result <- se_counts_weitrix(...)
    saveHDF5SummarizedExperiment(result, out_path, replace=TRUE)
    invisible()
})



do_peaks_shift <- function(
        out_path, peaks_se, organsism,
        mispriming=FALSE,utr_only=TRUE) {
}


#' @export
do_weitrix_components <- function(out_path, weitrix, p=30, ...) working_in(out_path, {
    weitrix <- load_banquet(weitrix)

    comp_seq <- weitrix_components_seq(weitrix, p=p, ...)
    saveRDS(comp_seq, file.path(out_path,"comp_seq.rds"))

    rand_weitrix <- weitrix_randomize(weitrix)
    rand_comp_seq <- weitrix_components_seq(rand_weitrix, p=p, ...)
    saveRDS(rand_comp_seq, file.path(out_path,"rand_comp_seq.rds"))

    invisible()
})






