#
# Automated analysis steps, producing banquetry
#
# do_ functions create a file or directory rather than returning a value.
# Their inputs can also be the path to a file or directory.
#


#' @export
PIPELINE_STAGES = c(
    "load",
    "assign",
    "weitrices",
    "report"
)


#' Analyse the output of polyApipe.py
#'
#' Analyse the output of polyApipe.py. The resulting directory can be loaded
#' with \code{load_banquet}. An HTML report is also produced.
#'
#' You can choose to run only specific stages of the pipeline with the \code{stages} argument. See the global \code{PIPELINE_STAGES} for possible stages. \code{stages} is a character vector or an integer vector indexing into \code{PIPELINE_STAGES}.
#'
#' @param out_path Output directory name.
#'
#' @param counts_file_dir A directory containing counts files. Batch names are taken from the basename of the count files without the .tab.gz suffix. Give either this argument or \code{counts_files} but not both.
#'
#' @param counts_files Alternative to counts_file_dir. One or more filenames for .tab.gz files produced by polyApipe.py. Give either this argument or \code{counts_file_dir} but not both.
#'
#' @param batch_names If using counts_files, give a vector of the batch/sample names, in same order as counts_files
#'
#' @param peak_info_file GTF formatted peak file as output from polyApipe.py.
#'
#' @param organism Organism directory, as created by \code{do_ensembl_organism}.
#'
#' @param remove_mispriming Remove peaks considered to be mispriming peaks.
#'
#' @param utr_or_extension_only Remove peaks in the 5'UTR, exons, or introns.
#'
#' @param cells_to_use Character vector of cells to use.
#'
#' @param peak_min_present A peak is retained if it is present in this number of cells.
#' 
#' @param peak_min_prop A peak is used for APA and gene expression calculations only if it constitutes at least this proportion of total UMIs for the gene.
#'
#' @param do_logNormCounts Compute logcounts for gene and peak expression?
#'
#' @param do_computeSumFactors When computing logcounts, use scran::computeSumFactors? If not, unadjusted cell library sizes are used.
#'
#' @param title Title to use in report.
#'
#' @param stages Stages of the pipeline to run. Stages are listed in PIPELINE_STAGES.
#'
#' @return 
#' There is no return value. Results are placed in the \code{out_path} directory. They can be loaded using \code{load_banquet(out_path)}.
#'
#' @export
do_pipeline <- function(
        out_path,
        
        # Either:
        counts_file_dir=NULL,
        # Or:
        counts_files=NULL,
        batch_names="",
        
        peak_info_file,
        organism,

        cell_name_func=function(batch, cell) paste0(batch, cell),
        cells_to_use=NULL,

        # peaks_weitrices options
        remove_mispriming=TRUE, 
        utr_or_extension_only=FALSE,
        peak_min_present=50,
        peak_min_prop=0.01,
        # log expression level options
        do_logNormCounts=TRUE,
        do_computeSumFactors=TRUE, 
                
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
        message("-- 1/4 load --")
    
        if (!is.null(counts_file_dir)) {
            assert_that(is.null(counts_files))
            
            load_peaks_counts_dir_into_sce(
                counts_file_dir=counts_file_dir,
                peak_info_file=peak_info_file,
                output=peak_counts,
                cell_name_func=cell_name_func,
                cells_to_use=cells_to_use,
                replace=TRUE)
        } else {
            load_peaks_counts_files_into_sce(
                counts_file_list=counts_files,
                batch_names=batch_names,
                peak_info_file=peak_info_file,
                output=peak_counts,
                cell_name_func=cell_name_func,
                cells_to_use=cells_to_use,
                replace=TRUE)
        }
    }

    if ("assign" %in% stages) {
        message("-- 2/4 assign --")
        do_se_reassign(peak_counts, organism)
        
        write_gff3(
            rowRanges(load_banquet(peak_counts)),
            file.path(out_path, "peaks.gff3"), 
            index=TRUE)
    }
    
    if ("weitrices" %in% stages) {
        message("-- 3/4 weitrices --")
        do_peaks_weitrices(out_path, peak_counts,
            peak_min_present=peak_min_present,
            peak_min_prop=peak_min_prop,
            cells_to_use=cells_to_use,
            remove_mispriming=remove_mispriming, 
            utr_or_extension_only=utr_or_extension_only,
            do_logNormCounts=do_logNormCounts,
            do_computeSumFactors=do_computeSumFactors)
    }
    
    if ("report" %in% stages) {
        message("-- 4/4 report --")
        do_pipeline_report(
            out_path,
            title=paste(title,"- overview"))
    }

    invisible(NULL)
}



# Add logcounts and some summary information to a SingleCellExperiment object
se_logcounts <- function(se, size_factors=NULL, do_logNormCounts=TRUE, do_computeSumFactors=TRUE) {
    se <- load_banquet(se)
    
    col_sums <- colSums(assay(se,"counts"))
    row_sums <- rowSums(assay(se,"counts"))

    col_present <- colSums(assay(se,"counts") > 0)
    row_present <- rowSums(assay(se,"counts") > 0)

    if (do_logNormCounts) {
        if (do_computeSumFactors && is.null(size_factors)) {
            message("Preclustering")
            preclusters <- quickCluster(se, block=colData(se)$batch, block.BPPARAM=bpparam())
            message(length(unique(preclusters)), " preclusters")
            colData(se)$preclusters <- preclusters
            se <- computeSumFactors(se, cluster=preclusters, BPPARAM=bpparam())
        } else {
            if (is.null(size_factors))
                size_factors <- col_sums 
            sizeFactors(se) <- size_factors
        }

        # Prevent errors for any completely missing cells
        sizeFactors(se)[col_present==0] <- 1
    
        se <- logNormCounts(se)
    }

    rowData(se)$total_count <- row_sums
    rowData(se)$total_present <- row_present
    colData(se)$total_count <- col_sums
    colData(se)$total_present <- col_present
    
    se
}


amalgamate <- function(vec) paste(unique(vec),collapse="/")

realize_RleArray <- function(mat) realize(mat, BACKEND="RleArray")

#' Shifts and proportions weitrices from peak counts
#'
#' Convert a SummarizedExperiment of peak counts into a weitrix of shifts and a weitrix of proportions.
#'
#' Peaks with less than an average of min_peak_count_avg counts will be discarded.
#'
#' Only genes with two or more peaks are included.
#'
#' @export
do_peaks_weitrices <- function(
        out_path,
        peaks_se,
        remove_mispriming=TRUE, 
        utr_or_extension_only=FALSE,
        cells_to_use=NULL,
        peak_min_present=50,
        peak_min_prop=0.01,
        # log expression level options
        do_logNormCounts=TRUE,
        do_computeSumFactors=TRUE) {

    peaks_se <- load_banquet(peaks_se)

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

    # Remove unwanted cells and low count peaks
    good_cells <- cells_to_use[cells_to_use %in% colnames(peaks_se_relevant)]
    assert_that(length(good_cells) >= 1)
    
    peak_presence <- rowSums(
        assay(peaks_se_relevant,"counts")[,good_cells,drop=F] 
        > 0)
    good_peaks <- peak_presence >= peak_min_present
    
    peaks_se_final <- peaks_se_relevant[good_peaks, good_cells]
    # May be faster, but potentially uses lots of memory.
    #assay(peaks_se_final,"counts") <- realize(assay(peaks_se_final,"counts"))
    rowData(peaks_se_final)$total_count <- rowSums(assay(peaks_se_final,"counts"))

    message("Kept ", nrow(peaks_se_final), " peaks in ", ncol(peaks_se_final), " cells")

    # Work out peak grouping

    rd <- rowRanges(peaks_se_final) %>%
        as.data.frame() %>%
        rownames_to_column("name") %>%
        arrange(symbol_unique, seqnames, strand, ifelse(strand == "-", -start, end))
    # Grouped symbol_unique
    # Leftovers grouped by chromosome, strand
    # Within this, by position on strand

    # Apply min_prop filter
    rd_with_gene <- rd %>%
        filter(!is.na(symbol_unique)) %>%
        group_by(symbol_unique) %>%
        filter( total_count >= sum(total_count)*peak_min_prop ) %>%
        ungroup()
    
    rowData(peaks_se_final)$passed_min_prop <- rownames(peaks_se_final) %in% rd_with_gene$name

    # Genes to use for APA shift    
    rd2 <- rd_with_gene %>%
        group_by(symbol_unique) %>%
        filter(length(name) >= 2) %>%
        ungroup()

    gene_info <- rd_with_gene %>%
        group_by(symbol_unique) %>%
        summarize(
            gene_id=amalgamate(gene_id),
            symbol=amalgamate(symbol),
            biotype=amalgamate(biotype),
            peaks=paste(name,collapse=":"),
            regions=paste(region_grouped,collapse=":"))

    grouping <- rd2 %>%
        select(group=symbol_unique, name=name)

    message(nrow(gene_info), " genes")
    message(length(unique(grouping$group)), " genes with two or more peaks")
    
    # Outputs based on one or more peaks
    
    gene_counts <- rowsum(
        assay(peaks_se_final,"counts")[rd_with_gene$name,,drop=F], 
        rd_with_gene$symbol_unique)
    
    result_gene <- SingleCellExperiment(list(counts=gene_counts))
    colData(result_gene) <- colData(peaks_se_final)
    rowData(result_gene) <- 
        gene_info[match(rownames(result_gene), gene_info$symbol_unique),,drop=F]
    
    message("Compute normalized, log transformed counts")
    result_gene <- se_logcounts(result_gene, 
        do_logNormCounts=do_logNormCounts, do_computeSumFactors=do_computeSumFactors)

    saveHDF5SummarizedExperiment(
        result_gene, file.path(out_path,"gene_expr"), replace=TRUE)
    gene_size_factors <- sizeFactors(result_gene)
    rm(result_gene)

    # Use rd ordering    
    result_peak <- se_logcounts(peaks_se_final[rd$name,],
        size_factors=gene_size_factors, do_logNormCounts=do_logNormCounts, do_computeSumFactors=FALSE)

    saveHDF5SummarizedExperiment(
        result_peak, file.path(out_path,"peak_expr"), replace=TRUE)
    rm(result_peak)


    # Outputs based on two or more peaks
            
    message("Compute APA shifts and peak proportions")
    result_shift <- counts_shift(assay(peaks_se_final,"counts"), grouping, typecast=realize_RleArray)

    gene_info <-
        gene_info[match(rownames(result_shift), gene_info$symbol_unique),,drop=F]
    rowData(result_shift)$gene_id <- gene_info$gene_id
    rowData(result_shift)$symbol <- gene_info$symbol
    rowData(result_shift)$regions <- gene_info$regions
    rowData(result_shift)$biotype <- gene_info$biotype
    colData(result_shift) <- colData(peaks_se_final)

    saveHDF5SummarizedExperiment(
        result_shift, file.path(out_path,"shift"), replace=TRUE)
    rm(result_shift)

    result_prop <- counts_proportions(assay(peaks_se_final,"counts"), grouping, typecast=realize_RleArray)
    rowData(result_prop) <- rowData(peaks_se_final)[rownames(result_prop),]
    colData(result_prop) <- colData(peaks_se_final)

    saveHDF5SummarizedExperiment(
        result_prop, file.path(out_path,"prop"), replace=TRUE)
    rm(result_prop)

    out_path
}






