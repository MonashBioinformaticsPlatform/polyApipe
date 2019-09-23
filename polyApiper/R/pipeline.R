#
# Automated analysis steps, producing banquetry
#
# do_ functions create a file or directory rather than returning a value.
# Their inputs can also be the path to a file or directory.
#

#' @rdname do_pipeline
#' @export
PIPELINE_STAGES = c(
    "load",
    "assign",
    "weitrices",
    "shift_comp",
    "gene_expr_comp",
    "report"
)


#' Analyse the output of polyApipe.py
#'
#' Analyse the output of polyApipe.py. The resulting directory can be loaded
#' with \code{load_banquet}. An HTML report is also produced.
#'
#' You can choose to run only specific stages of the pipeline with the \code{stages} argument, which is a character vector. See the global \code{PIPELINE_STAGES} for possible stages.
#'
#' @param out_path Output directory name.
#'
#' @param counts_files One or more filenames for .tab.gz files produced by polyApipe.py. Give either this argument or \code{counts_file_dir} but not both.
#'
#' @param batch_names A vector of the batch/sample names, in same order as counts_files
#'
#' @param counts_file_dir A directory containing counts files. Give either this argument or \code{counts_files} but not both.
#'
#' @param peak_info_file GTF formatted peak file as output from polyApipe.py.
#'
#' @param organism Organism directory, as created by \code{do_ensembl_organism}.
#'
#' @param downsample If less than 1, downsample cells to this proportion when creating weitrices.
#'
#' @param remove_mispriming Remove peaks considered to be mispriming peaks.
#'
#' @param utr_or_extension_only Remove peaks in the 5'UTR, exons, or introns.
#'
#' @param present_min_count A gene is considered present in a cell with this count.
#'
#' @param min_cell_present_vs_avg Controls filtering of cells. Cells are retained if this proportion of the average number of present peaks are present. This is iterated, so the average is over the retained cells.
#'
#' @param min_peak_present After filtering cells, a peak is retained if it is present in this proportion of cells.
#' 
#' @param do_computeSumFactors Use scran::computeSumFactors? If not, unadjusted cell library sizes are used.
#'
#' @param p Number of components to find.
#'
#' @param n_restarts Use this many restarts when finding components. If the scree plot is uneven, increasing this may produce better results.
#'
#' @param title Title to use in report.
#'
#' @param stages Stages of the pipeline to run. Stages are listed in PIPELINE_STAGES.
#'
#' @return 
#' There is no return value. Results are placed in the \code{out_path} directory.
#'
#' @export
do_pipeline <- function(
        out_path,
        
        # Either:
        counts_files=NULL,
        batch_names="",
        # Or:
        counts_file_dir=NULL,
        
        peak_info_file,
        organism,
        
        # peaks_weitrices options
        downsample=1,
        remove_mispriming=TRUE, 
        utr_or_extension_only=FALSE,
        present_min_count=1,
        min_cell_present_vs_avg=0.5,
        min_peak_present=0.01,
        # log expression level options
        do_computeSumFactors=TRUE, 
        
        # finding components
        p=20,
        n_restarts=2,
        
        title="polyApiper pipeline run",
                
        stages=1:length(PIPELINE_STAGES)) {

    if (is.numeric(stages)) 
        stages <- PIPELINE_STAGES[stages]
    assert_that(all(stages %in% PIPELINE_STAGES))

    ensure_dir(out_path)

    # Record invocation
    env <- environment() 
    cat("do_pipeline(\n",
        map_chr(names(formals()), 
            ~paste0("    ",.,"=",
                paste(deparse(env[[.]]), collapse="\n        ") )) %>%
            paste(collapse=",\n"),
        ")\n",
        sep="",file=file.path(out_path,"invocation.txt"))


    peak_counts <- file.path(out_path, "raw_peak_counts")

    if ("load" %in% stages) {
        message("-- load --")
    
        if (!is.null(counts_file_dir)) {
            assert_that(is.null(counts_files))
            
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
            downsample=downsample,
            present_min_count=present_min_count,
            min_peak_present=min_peak_present,
            min_cell_present_vs_avg=min_cell_present_vs_avg,
            remove_mispriming=remove_mispriming, 
            utr_or_extension_only=utr_or_extension_only,
            do_computeSumFactors=do_computeSumFactors)
    }

    if ("shift_comp" %in% stages) {
        message("-- shift_components --")
        do_weitrix_components(
            file.path(out_path, "shift"),
            file.path(out_path, "shift","weitrix"),
            design=~0, p=p, n_restarts=n_restarts)
    }
    
    if ("gene_expr_comp" %in% stages) {
        message("-- expr_components --")
        do_weitrix_components(
            file.path(out_path, "gene_expr"),
            file.path(out_path, "gene_expr","weitrix"),
            p=p, n_restarts=n_restarts)
    }
    
    if ("report" %in% stages) {
        message("-- report --")
        do_pipeline_report(
            out_path,
            title=paste(title,"- overview"))
    }

    invisible(NULL)
}




#' @export
se_counts_weitrix <- function(se, do_computeSumFactors=TRUE) {
    se <- load_banquet(se)
    
    if (do_computeSumFactors) {
        se <- computeSumFactors(se) 
        #Could specify this, but it crashed when I tested it:
        # BPPARAM=getAutoBPPARAM()
    } else
        sizeFactors(se) <- colSums(assay(se,"counts"))
    
    se <- normalize(se)

    assay(se, "weights") <- matrix(1,nrow=nrow(se),ncol=ncol(se))
    result <- bless_weitrix(se, "logcounts", "weights")
    
    rowData(result)$total_reads <- rowSums(assay(se,"counts"))
    metadata(result)$weitrix$trend_formula <- "~splines::ns(log(total_reads),3)"
    
    result
}


#' @export
do_se_counts_weitrix <- function(out_path, ...) working_in(out_path, {
    result <- se_counts_weitrix(...)

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
        downsample=1,
        remove_mispriming=TRUE, 
        utr_or_extension_only=FALSE,
        present_min_count=1,
        min_cell_present_vs_avg=0.5,
        min_peak_present=0.05,
        # log expression level options
        do_computeSumFactors=TRUE) {
    peaks_se <- load_banquet(peaks_se)
    
    if (downsample < 1) {
       message("Downsampling ", downsample, " from ", ncol(peaks_se), " cells")
       subsample <- runif(ncol(peaks_se)) <= downsample
       peaks_se <- peaks_se[, subsample]
    }

    message(nrow(peaks_se), " peaks in ", ncol(peaks_se), " cells")

    # Keep relevant peaks
    
    keep <- rep(TRUE, nrow(peaks_se))    
    if (remove_mispriming)
        keep <- keep & !rowData(peaks_se)$misprime
    if (utr_or_extension_only)
        keep <- keep & (
            rowData(peaks_se)$region == "3'UTR" |
            rowData(peaks_se)$region == "extension")

    peaks_se_relevant <- peaks_se[keep,]

    # Remove low count cells and peaks

    cell_present <- colSums(assay(peaks_se_relevant,"counts") >= present_min_count)
    
    # Iteratively trim to a good set of cells
    good_cells <- rep(TRUE,length(cell_present))
    this_n_good <- sum(good_cells)
    repeat {
        min_cell <- max(1, mean(cell_present[good_cells]) * min_cell_present_vs_avg)
        good_cells <- cell_present >= min_cell
        last_n_good <- this_n_good
        this_n_good <- sum(good_cells)
        if (this_n_good == last_n_good) break;
    } 
    
    peak_presence <- rowMeans(
        assay(peaks_se_relevant,"counts")[,good_cells,drop=F] 
        >= present_min_count)
    good_peaks <- peak_presence >= min_peak_present
    
    peaks_se_final <- peaks_se_relevant[good_peaks, good_cells]
    assay(peaks_se_final,"counts") <- realize(assay(peaks_se_final,"counts"))

    message("Kept ", nrow(peaks_se_final), " peaks in ", ncol(peaks_se_final), " cells")

    # Work out peak grouping

    rd <- rowRanges(peaks_se_final) %>%
        as.data.frame() %>%
        rownames_to_column("name") %>%
        arrange(gene_id, seqnames, strand, ifelse(strand == "-", -start, end))
    # Grouped gene_id
    # Leftovers grouped by chromosome, strand
    # Within this, by position on strand
    
    rd_with_gene <- rd %>%
        filter(!is.na(gene_id))
    
    rd2 <- rd_with_gene %>%
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
    
    gene_counts <- realize(rowsum(
        assay(peaks_se_final,"counts")[rd_with_gene$name,,drop=F], 
        rd_with_gene$gene_id))
    
    result_gene <- SingleCellExperiment(list(counts=gene_counts))
    colData(result_gene) <- colData(peaks_se_final)
    rowData(result_gene) <- 
        gene_info[match(rownames(result_gene), gene_info$gene_id),,drop=F]
    
    message("Compute normalized, log transformed counts")
    result_gene <- se_counts_weitrix(result_gene, 
        do_computeSumFactors=do_computeSumFactors)
    
    # Use rd ordering    
    result_peak <- se_counts_weitrix(peaks_se_final[rd$name,],
        do_computeSumFactors=do_computeSumFactors)
    
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






