#' Load peak counts file into a SingleCellExperiment object
#'
#' Load one 3-column counts.tab.gz file from polyApipe.py into a 
#' SingleCellExperiment for analysis.
#'
#' @param counts_file Tab-delimeted file of read counts per peak. 
#' The '<sample>.tab.gz' as output from polyApipe.py. 
#' Format: <peak> <cell> <count>
#' @param peak_info_file GTF formatted peak file _specfically_ as output from polyApipe.py
#' @param output Where to save the sce object on disk using ssaveHDF5SummarizedExperiment.
#' Should not exist. This directory will be created.
#' @param batch A name of the batch/sample. Important if multiple samples are 
#' going to be combined.
#' @param cell_name_func A function taking two arguments, the batch name and a vector of cell barcodes. Returns a vector of cell names. Defaults to just paste0-ing the batch and cell barcode together.
#' @param cells_to_use A vector of cell names to use (should match the output of cell_name_func). NULL to retain all cells.
#' @param n_max Only read the first *n_max* lines of counts table. (Default = -1 (all))
#' @param ...  Extra options passed to saveHDF5SummarizedExperiment 
#' 
#' @examples  
#' 
#' counts_file    <- system.file("extdata", "demo_dataset/demo_counts/SRR5259354_demo.tab.gz", package = "polyApiper")
#' peak_info_file <- system.file("extdata", "demo_dataset/demo_polyA_peaks.gff", package = "polyApiper") 
#' 
#' sce <- load_peaks_counts_into_sce(counts_file    =counts_file1, 
#'                                   peak_info_file = peak_info_file, 
#'                                   output         = "demo") 
#' 
#' \dontrun{
#' 
#' sce <- load_peaks_counts_into_sce(real_counts_file, "real_polyA_peaks.gff", output = "mysce")
#' 
#' sce <- load_peaks_counts_into_sce(counts_file    = "run1_knockout_counts.tab.gz", 
#'                                   peak_info_file = "run1_knockout_polyA_peaks.gff", 
#'                                   output         = "./knockout_rep1", 
#'                                   batch          = "knockout_rep1")
#' 
#' # To reload from disk later 
#' library('HDF5Array')
#' library('SingleCellExperiment')
#' sce <- loadHDF5SummarizedExperiment("knockout_rep1/")
#' 
#' }
#' 
#' @family peak counts loading functions
#' 
#' @importFrom S4Vectors  DataFrame Rle 
#' @importFrom dplyr select
#' @importFrom utils read.table
#' @importFrom Matrix sparseMatrix Matrix
#' @importFrom DelayedMatrixStats colSums2
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom HDF5Array saveHDF5SummarizedExperiment writeTENxMatrix
#' @importFrom IRanges IRanges
#' @import GenomicRanges 
#' 
#' @export
load_peaks_counts_into_sce <- function(
      counts_file, 
      peak_info_file, 
      batch="",
      output,
      cell_name_func = function(batch, cell) paste0(batch, cell),
      cells_to_use = NULL,
      n_max=-1,
      missing_peaks_are_zero = TRUE,
      ... ) {
   
   message("Loading ", counts_file)

   if (output == "" | output == "." | output == "./") {
      stop("Specify some output name to write the sce hdf data, with ",
           "'output' parameter. A directory will be created there.")
   }
   
   # Get list of peaks from the peaks .gtf info file. 
   # Because no all files/cells will have all peaks with multipl samples.
   # need to have a consistant list anything not seen is assumed 0.
   # Obviously needs to be same file for files to merge, so order is the same.
   peak_info <- S4Vectors::DataFrame(
      read_polyA_peak_file_gtf(peak_info_file) %>%
      dplyr::select(peak, chr, start, end, strand, peakdepth, misprime)
      ) 
   rownames(peak_info) <- peak_info$peak
   peak_info$peakdepth <- as.numeric(peak_info$peakdepth)
   peak_info$misprime  <- as.logical((peak_info$misprime))
   
   #Then, store in granges object
   gr <- GRanges(
      seqnames = Rle(factor(peak_info$chr)),
      ranges   = IRanges(start=peak_info$start, end= peak_info$end),
      strand   = Rle(factor(peak_info$strand))
   )
   mcols(gr) <- peak_info[, ! colnames(peak_info) %in% c("chr", "start", "end", "strand")]
   names(gr) <- gr$peak
   
   
   # Not using read_tsv because of compression.(segfault.) 
   # NB: https://github.com/tidyverse/readr/issues/610
   the_counts <- read.table(counts_file, sep="\t", header=TRUE, 
                            as.is=TRUE, 
                            colClasses = c("factor","factor", "numeric"),
                            nrows=n_max)
   colnames(the_counts) <- c("peak","cell","count")
   
   # cell_name_func to make cell barcodes unique across multi-sample experiments.
   cell_to_barcode        <- as.character(the_counts$cell)
   names(cell_to_barcode) <- cell_name_func(batch, the_counts$cell)
   the_counts$cell        <- factor(names(cell_to_barcode))

   # Optionally only retain desired cells   
   if (!is.null(cells_to_use)) {
       the_counts <- the_counts[the_counts$cell %in% cells_to_use,]
       the_counts$cell <- droplevels(the_counts$cell)
   }
   
   # First construct the sparse matrix of all the counts
   counts_matrix <- Matrix::sparseMatrix(
      i = as.integer(the_counts$peak),
      j = as.integer(the_counts$cell),
      x = as.integer(the_counts$count))
   colnames(counts_matrix) <- levels(the_counts$cell)
   rownames(counts_matrix) <- levels(the_counts$peak)   

   message("Loaded ", nrow(counts_matrix), " x ", ncol(counts_matrix), " matrix of counts")
   
   # Assume that any gene not present has counts of 0, and match order.
   # (This is used to make consistant gene-list ses to merge later)
   full_peak_list  <- peak_info$peak 
   absent_genes    <- full_peak_list[ ! full_peak_list %in% rownames(counts_matrix) ]
   stopifnot( all(rownames(counts_matrix) %in% full_peak_list))
   
   counts_matrix.absent <- Matrix::Matrix(0, nrow=length(absent_genes), ncol=ncol(counts_matrix), sparse = TRUE)  
   rownames(counts_matrix.absent) <- absent_genes
   
   counts_matrix <- rbind(counts_matrix, counts_matrix.absent)
   counts_matrix <- counts_matrix[full_peak_list,]
   
   # cell data
   new_colData <- DataFrame(
      cell        = colnames(counts_matrix), 
      barcode     = cell_to_barcode[colnames(counts_matrix)], 
      batch       = batch) 
   
   # Load up the sce.
   sce <- SingleCellExperiment(
      assays    = list(counts = counts_matrix ),
      colData   = new_colData,
      rowRanges = gr[rownames(counts_matrix),] )
   rm(counts_matrix)
      
   # Save it
   # - As at Bioconductor 3.12, cbinding when loading multiple matrices 
   #   does not work unless as.sparse=FALSE given here.
   sce <- saveHDF5SummarizedExperiment(sce, 
                                       output, 
                                       as.sparse=FALSE,
                                       ... ) 
  
   return(sce)
}


#' Load multiple peak counts files into a SingleCellExperiment object
#'
#' Load multiple named 3-column counts.tab.gz files from polyApipe.py into a 
#' SingleCellExperiment for analysis. Each cell will be uniquely named, 
#' and colData will include the batch.
#' 
#' @param counts_file_list A vector of tab-delimeted files of read counts per peak. 
#' The '<sample>.tab.gz' as output from polyApipe.py. 
#' Format: <peak> <cell> <count>
#' @param peak_info_file GTF formatted peak file _specfically_ as output from polyApipe.py
#' @param output Where to save the sce object on disk using saveHDF5SummarizedExperiment.
#' Should not exist. This directory will be created. 
#' @param batch_names A vector of the batch/sample names, in same order as counts_file_list
#' @param replace  Overwrite existing output file? (Default=FALSE)
#' @param ... Other parameters passed through to *load_peaks_counts_into_sce* 
#' e.g. n_max, min_reads_per_barcode (but not cell_prefix or missing_peaks_are_zero)
#' 
#' 
#' @examples 
#' 
#' counts_file1   <- system.file("extdata", "demo_dataset/demo_counts/SRR5259354_demo.tab.gz", package = "polyApiper")
#' counts_file2   <- system.file("extdata", "demo_dataset/demo_counts/SRR5259422_demo.tab.gz", package = "polyApiper")
#' peak_info_file <- system.file("extdata", "demo_dataset/demo_polyA_peaks.gff", package = "polyApiper") 
#' 
#' sce <- load_peaks_counts_files_into_sce(c(counts_file1,counts_file2),
#'                batch_names    = c("demo1","demo2"),
#'                peak_info_file = peak_info_file,
#'                output = "demo1and2")
#' 
#' @family peak counts loading functions
#' 
#' @export
load_peaks_counts_files_into_sce <- function(
         counts_file_list, batch_names, # Order must match!
         peak_info_file, 
         output,
         replace=FALSE,
         ... ) {
   
   # Check names
   #NB: here, batchnames Will be used for cell_prefix nyah nyah 
   if (length(counts_file_list) != length(batch_names) | 
       length(counts_file_list) != length(unique(counts_file_list)) | 
       length(batch_names) != length(unique(batch_names)) ) {
          stop(c("File list and batch names list should be the same length, and without duplicates. ",
                 "Have batches: ", paste(batch_names,collapse=", "), 
                 ' for ',  paste(counts_file_list, collapse=", ")))
   }

   # Handle large datasets by running separately, and then joining.
   # interium ses need to be saved
   outputs <- paste0(tempfile(),"_",batch_names)
   sce <- do.call(cbind,
       mapply(FUN=load_peaks_counts_into_sce, 
                       counts_file   = counts_file_list,
                       batch         = batch_names,   
                       output = outputs, 
                       MoreArgs = list(
                          peak_info_file=peak_info_file,
                          missing_peaks_are_zero = TRUE,
                          #replace is not passed through, should be in tmp
                          ...
      )))

   # Save it (and use replace here.)
   sce <- saveHDF5SummarizedExperiment(sce, output, replace=replace)
   
   # Delete the temp files
   file.remove(file.path(outputs,"assays.h5"))
   file.remove(file.path(outputs,"se.rds"))
   unlink(outputs)
   
   return(sce)
}





#' Load a directory of peak counts files into a SingleCellExperiment object
#'
#' Load a directory's worth of 3-column counts.tab.gz files from polyApipe.py 
#' into a SingleCellExperiment for analysis. Each cell will be uniquely named, 
#' and colData will include the batch. The 'batch_regex' can be used to extract 
#' nice sample names from filenames.
#'  
#'
#' @param counts_file_dir The directory full of counts files (in .tab.gz) 
#' as output from polyApipe.py script (usually called <something>_counts)
#' @param peak_info_file GTF formatted peak file _specfically_ as output from polyApipe.py
#' @param output Where to save the sce object on disk using saveHDF5SummarizedExperiment.
#' Should not exist. This directory will be created. 
#' @param batch_regex A string of a regex to extract a nice sample name from a the counts 
#' files in counts_file_dir . 
#' This will become the batch name (and prefix of cell name) for each cell in colData.
#' If not specified this will just use the filename up to but excludeing '.tab.gz'.
#' The regex should have one capture group '()' and pull out something unique for every file.
#' (Default = '(.*)\\.tab\\.gz')
#' @param ... Other parameters passed through to 
#' load_peaks_counts_files_into_sce/load_peaks_counts_into_sce
#' 
#' @examples 
#' 
#' counts_dir     <- dirname(system.file("extdata", "demo_dataset/demo_counts/SRR5259354_demo.tab.gz", package = "polyApiper"))
#' peak_info_file <- system.file("extdata", "demo_dataset/demo_polyA_peaks.gff", package = "polyApiper") 
#' 
#' sce <- load_peaks_counts_dir_into_sce(counts_dir, 
#'                                peak_info_file = peak_info_file, 
#'                                output = "demodirsce",
#'                                min_reads_per_barcode=1) # Required for tiny demo data.
#'
#'
#' \dontrun{
#'
#' sce <- load_peaks_counts_dir_into_sce("myexpr_counts/", "my_expr_polyA_peaks.gff", "./myexpr_sce") 
#'
#' sce <- load_peaks_counts_dir_into_sce("data/MyExperimentFullRun20190516_counts/", 
#'                                peak_info_file = "data/MyExperimentFullRun20190516_polyA_peaks.gff", 
#'                                output         = "data/MyExperimentFullRun20190516_sce",
#'                                batch_regex    = 'FullRun_(.*)\\.tab\\.gz'
#'                                )
#' }
#' 
#'
#' @family peak counts loading functions
#'
#' @export
load_peaks_counts_dir_into_sce <- function(
      counts_file_dir, 
      peak_info_file,
      output,
      batch_regex = '(.*)\\.tab\\.gz',
      ... ) {
   
   counts_file_list <- list.files(counts_file_dir)
   batch_names      <-  sub( batch_regex,'\\1', basename(counts_file_list))
   
   sce <- load_peaks_counts_files_into_sce(
      counts_file_list   = file.path(counts_file_dir,counts_file_list), 
      batch_names        = batch_names, # Order must match!
      peak_info_file     = peak_info_file, 
      output    = output,
      ... )
   
   return(sce)
}



   


#' Read in the polyA peaks gtf file
#' 
#' The polyApipe.py script will write a gtf formatted file with the polyA peak 
#' locations. This function will read it into a tabular format in R. 
#' This is not a general gff/gtf reader funtion -  
#' it will break on probably any other gff!
#' 
#' @param peak_info_file GTF formatted peak file _specfically_ as output from polyApipe.py
#' 
#' @examples 
#' 
#' peak_info_file <- system.file("extdata", "demo_dataset/demo_polyA_peaks.gff", package = "polyApiper") 
#' peak_info_table <- read_polyA_peak_file_gtf(peak_info_file)
#' 
#' 
#' @importFrom readr read_tsv
#' @importFrom stringr str_trim
#' @importFrom tibble tibble
#' @importFrom tidyr separate_rows separate spread
#' 
#' @export
read_polyA_peak_file_gtf <- function(gtf_file) {
   # Very specifically this format of gtf.
   #Y	polyAends	polyApeak	20592383	20592633	.	+	.	peak="Y_20592632"_f; peakdepth="64"; misprime="False";
   #1	polyAends	polyAends	3630788	3631037	.	-	.  peak="1_3630788_r"; peakdepth="256"; misprime="False";
   raw_gtf_table <- read_tsv(gtf_file, comment = "#", col_types = 'cccii?ccc', na=".",
                             col_names = c("chr","source","feature","start","end","score","strand","frame","attribute" )) 
   raw_gtf_table$attribute <- gsub('"','', str_trim(raw_gtf_table$attribute))
   
   #peak="1_3630788_r"; peakdepth="256"; misprime="False";
   att_table <- tibble(attribute = raw_gtf_table$attribute, splitter=raw_gtf_table$attribute) %>% separate_rows(sep='; ', splitter)
   att_table$splitter <- gsub(';$', '', str_trim(att_table$splitter) )
   att_table <- dplyr::filter(att_table , splitter != "")
   
   att_table <- att_table %>% 
      separate(splitter,c("att_name","att_value"), sep="=") %>%
      spread(key=att_name, value=att_value)
   
   # Now attributes are columns
   polyA_peak_feature_table <- merge(raw_gtf_table, att_table, all.x = TRUE, by="attribute" )
   
   stopifnot(nrow(polyA_peak_feature_table)==nrow(raw_gtf_table), 
             nrow(polyA_peak_feature_table)==nrow(att_table))
   
   
   return(polyA_peak_feature_table[,-1])
}







