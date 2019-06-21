#
# Automated analysis steps, producing banquetry
#
# do_ functions create a file or directory rather than returning a value.
# Their inputs can also be the path to a file or directory.
#

#' @rdname do_pipeline
#' @export
pipeline_stages = c(
    "load",
    "assign",
    "weitrices",
    "gene_expr_comp",
    "shift_comp"
)

#' Analyse the output of polyApipe.py
#'
#' Analyse the output of polyApipe.py. The resulting directory can be loaded
#' with load_banquet.
#'
#' @export
do_pipeline <- function(
        out_path,
        
        # Either:
        counts_files,
        batch_names="",
        # Or:
        counts_file_dir,
        
        peak_info_file,
        organism,
        
        # peaks_weitrices options
        min_peak_count_avg=0.01,
        min_cell_count_vs_avg=0.25,
        remove_mispriming=TRUE, 
        utr_or_extension_only=FALSE,
        # log expression level options
        do_computeSumFactors=FALSE, 
        do_vooma=FALSE,
        
        # finding components
        p=20,
                
        stages=pipeline_stages) {

    ensure_dir(out_path)

    peak_counts <- file.path(out_path, "raw_peak_counts")

    if ("load" %in% stages) {
        message("-- load --")
    
        if (!missing(counts_file_dir)) {
            assert_that(missing(counts_file_list))
            
            load_peaks_counts_dir_into_sce(
                counts_file_dir=counts_file_dir,
                peak_info_file=peak_info_file,
                output=peak_counts,
                replace=TRUE)
        } else {
            load_peaks_counts_files_into_sce(
                counts_file_list=counts_files,
                batch_names=batch_names,
                peak_info_file=peak_info_file,
                output=peak_counts,
                replace=TRUE)
        }
    }

    if ("assign" %in% stages) {
        message("-- assign --")
        do_se_reassign(peak_counts, organism)
    }
    
    if ("weitrices" %in% stages) {
        message("-- weitrices --")
        do_peaks_weitrices(out_path, peak_counts,
            min_peak_count_avg=min_peak_count_avg,
            min_cell_count_vs_avg=min_cell_count_vs_avg,
            remove_mispriming=remove_mispriming, 
            utr_or_extension_only=utr_or_extension_only,
            do_computeSumFactors=do_computeSumFactors, 
            do_vooma=do_vooma)
    }
    
    if ("gene_expr_comp" %in% stages) {
        message("-- expr_components --")
        do_weitrix_components(
            file.path(out_path, "gene_expr"),
            file.path(out_path, "gene_expr","weitrix"),
            p=p)
    }

    if ("shift_comp" %in% stages) {
        message("-- expr_components --")
        do_weitrix_components(
            file.path(out_path, "shift"),
            file.path(out_path, "shift","weitrix"),
            p=p)
    }

    invisible(NULL)
}




#' @export
se_counts_weitrix <- function(se, do_computeSumFactors=FALSE, do_vooma=FALSE, plot_filename=NULL) {
    se <- load_banquet(se)
    message(nrow(se), " x ", ncol(se), " counts")
    
    if (do_computeSumFactors)
        se <- computeSumFactors(se)
    else
        sizeFactors(se) <- colSums(assay(se,"counts"))
    
    se <- normalize(se)

    if (do_vooma) {
        if (!is.null(plot_filename)) png(plot_filename)
        # TODO: better design matrix
        # TODO: low-memory solution
        result <- as_weitrix(vooma(assay(se,"logcounts"), plot=TRUE, span=0.25))
        rowData(result) <- rowData(se)
        colData(result) <- colData(se)
        if (!is.null(plot_filename)) dev.off()
        return(result)
    } else {
        assay(se, "weights") <- matrix(1,nrow=nrow(se),ncol=ncol(se))
        return(bless_weitrix(se, "logcounts", "weights"))
    }
}

#' @export
do_se_counts_weitrix <- function(out_path, ...) working_in(out_path, {
    result <- se_counts_weitrix(plot_filename=file.path(out_path,"vooma.png"), ...)

    saveHDF5SummarizedExperiment(result, file.path(out_path,"weitrix"), replace=TRUE)
    invisible()
})


amalgamate <- function(vec) paste(unique(vec),collapse="/")

#' Shifts and proportions weitrices from peak counts
#'
#' Convert a SummarizedExperiment of peak counts into a weitrix of shifts and a weitrix of proportions.
#'
#' Cells with less than min_cell_count_vs_avg times the average cell count will be discarded.
#'
#' Peaks with less than an average of min_peak_count_avg counts will be discarded.
#'
#' Only genes with two or more peaks are included.
#'
#' @export
peaks_weitrices <- function(
        peaks_se,
        min_peak_count_avg=0.01,
        min_cell_count_vs_avg=0.25,
        remove_mispriming=TRUE, 
        utr_or_extension_only=FALSE,
        # log expression level options
        do_computeSumFactors=FALSE, 
        do_vooma=FALSE
        ) {
    peaks_se <- load_banquet(peaks_se)

    message(nrow(peaks_se), " peaks in ", ncol(peaks_se), " cells")

    # Keep relevant peaks
    
    keep <- rep(TRUE, nrow(peaks_se))    
    if (remove_mispriming)
        keep <- keep & rowData(peaks_se)$misprime
    if (utr_or_extension_only)
        keep <- keep & (
            rowData(peaks_se)$region == "3'UTR" |
            rowData(peaks_se)$region == "extension")

    peaks_se_relevant <- peaks_se[keep,]

    # Remove low count cells and peaks

    cell_counts <- colSums(assay(peaks_se_relevant,"counts"))
    min_cell <- max(1, mean(cell_counts) * min_cell_count_vs_avg)
    good_cells <- cell_counts >= min_cell
    
    peak_means <- rowMeans(assay(peaks_se_relevant,"counts")[,good_cells,drop=F])
    good_peaks <- peak_means >= min_peak_count_avg
    
    peaks_se_final <- peaks_se_relevant[good_peaks, good_cells]

    message("Kept ", nrow(peaks_se_final), " peaks in ", ncol(peaks_se_final), " cells")

    # Work out peak grouping

    rd <- rowRanges(peaks_se_final) %>%
        as.data.frame() %>%
        rownames_to_column("name") %>%
        arrange(ifelse(strand == "-", -start, end))
    
    rd2 <- rd %>%
        filter(!is.na(gene_id)) %>%
        group_by(gene_id) %>%
        filter(length(name) >= 2) %>%
        ungroup()

    gene_info <- rd %>%
        filter(!is.na(gene_id)) %>%
        group_by(gene_id) %>%
        summarize(
            symbol=amalgamate(symbol),
            biotype=amalgamate(biotype),
            peaks=paste(name,collapse=":"),
            regions=paste(region,collapse=":"))

    grouping <- rd2 %>%
        select(group=gene_id, name=name, region=region)

    message(nrow(gene_info), " genes")
    message(length(unique(grouping$group)), " genes with two or more peaks")
    
    # Outputs based on one or more peaks
    
    good <- !is.na(rd$gene_id)
    gene_counts <- rowsum(assay(peaks_se_final,"counts")[good,,drop=F], rd$gene_id[good])
    
    result_gene <- SingleCellExperiment(list(counts=gene_counts))
    colData(result_gene) <- colData(peaks_se_final)
    rowData(result_gene) <- 
        gene_info[match(rownames(result_gene), gene_info$gene_id),,drop=F]
    
    message("Compute normalized, log transformed counts")
    result_gene <- se_counts_weitrix(result_gene, 
        do_computeSumFactors=do_computeSumFactors, do_vooma=do_vooma)
    result_peak <- se_counts_weitrix(peaks_se_final,
        do_computeSumFactors=do_computeSumFactors, do_vooma=do_vooma)
    
    # Outputs based on two or more peaks
            
    message("Compute APA shifts and peak proportions")
    result_shift <- counts_shift(assay(peaks_se_final,"counts"), grouping)
    result_prop <- counts_proportions(assay(peaks_se_final,"counts"), grouping)

    gene_info <-
        gene_info[match(rownames(result_shift), gene_info$gene_id),,drop=F]
    rowData(result_shift)$symbol <- gene_info$symbol
    rowData(result_shift)$regions <- gene_info$regions
    rowData(result_shift)$biotype <- gene_info$biotype

    rowData(result_prop) <- rowData(peaks_se_final)[rownames(result_prop),]

    colData(result_shift) <- colData(peaks_se_final)
    colData(result_prop) <- colData(peaks_se_final)

    list(shift=result_shift, prop=result_prop, peak=result_peak, gene=result_gene)
}

#' @export
do_peaks_weitrices <- function(out_path, ...) working_in(out_path, {
    result <- peaks_weitrices(...)
    
    ensure_dir(file.path(out_path, "shift"))
    saveHDF5SummarizedExperiment(
        result$shift, file.path(out_path,"shift","weitrix"), replace=TRUE)

    ensure_dir(file.path(out_path, "prop"))
    saveHDF5SummarizedExperiment(
        result$prop, file.path(out_path,"prop","weitrix"), replace=TRUE)

    ensure_dir(file.path(out_path, "peak_expr"))
    saveHDF5SummarizedExperiment(
        result$peak, file.path(out_path,"peak_expr","weitrix"), replace=TRUE)

    ensure_dir(file.path(out_path, "gene_expr"))
    saveHDF5SummarizedExperiment(
        result$gene, file.path(out_path,"gene_expr","weitrix"), replace=TRUE)

    invisible()
})


#' @export
do_weitrix_components <- function(out_path, weitrix, p=20, ...) working_in(out_path, {
    weitrix <- load_banquet(weitrix)

    comp_seq <- weitrix_components_seq(weitrix, p=p, ...)
    saveRDS(comp_seq, file.path(out_path,"comp_seq.rds"))

    rand_weitrix <- weitrix_randomize(weitrix)
    rand_comp <- weitrix_components(rand_weitrix, p=1, ...)
    saveRDS(rand_comp, file.path(out_path,"rand_comp.rds"))

    invisible()
})






