
#' @importFrom readr read_tsv
#' @importFrom stringr str_trim
#' @importFrom tibble tibble
#' @importFrom tidyr separate_rows separate spread
#' 
#' @export
read_polyA_peak_file_gtf <- function(gtf_file, sep_genes=FALSE) {
   # Very specifically this format of gtf.
   #Y	polyAends	polyApeak	20592383	20592633	.	+	.	peakgene="EIF1AY"; peak="EIF1AY:20592632"; peakdepth="64"; misprime="False";
   #Y	polyAends	polyApeak	20843093	20843343	.	-	.	peakgene="Unknown"; peak="Unknown-Y-20843092"; peakdepth="14"; misprime="True";
   #Y	polyAends	polyApeak	21257670	21257920	.	+	.	peakgene="Unknown"; peak="Unknown-Y-21257669"; peakdepth="12"; misprime="False";
   #1	polyAends	polyAends	3630788	3631037	.	-	.	peakgene="Unknown"; peak="Unknown:1_3630788_r"; peakdepth="256"; misprime="False";
   raw_gtf_table <- read_tsv(gtf_file, comment = "#", col_types = 'cccii?ccc', na=".",
                             col_names = c("chr","source","feature","start","end","score","strand","frame","attribute" )) 
   raw_gtf_table$attribute <- gsub('"','', str_trim(raw_gtf_table$attribute))
   
   
   #peakgene="EIF1AY"; peak="EIF1AY:20592632"; peakdepth="64"; misprime="False";
   #peakgene=RP4-635E18.9; peak=RP4-635E18.8,RP4-635E18.9:11030527; peakdepth=16; misprime=False;
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
   
   
   # Repeat peaks when assigned to multiple genes?
   # (default no - so peakgene can be ; separated list of multiple genes.)
   if (sep_genes) {
      polyA_peak_feature_table <- separate_rows(polyA_peak_feature_table, peakgene, sep=",")
   }
   
   return(polyA_peak_feature_table[,-1])
}




#' I think ther was some reason to use DataFrame??
#' 
#' ...  Extra options passed to saveHDF5SummarizedExperiment E.g. verbose, replace=TRUE
#' 
#' 
#' @importFrom S4Vectors  DataFrame
#' @importFrom dplyr select
#' 
#' @importFrom utils read.table
#' @importFrom Matrix sparseMatrix Matrix
#' @importFrom DelayedMatrixStats colSums2
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom HDF5Array saveHDF5SummarizedExperiment writeTENxMatrix
#' 
#' 
#' @export
load_peaks_counts_into_sce.NOh5step <- function(counts_file, 
                                                peak_info_file, 
                                                batch="",
                                                sce_save_path = batch,
                                                cell_prefix   = batch,
                                                n_max=-1,
                                                min_reads_per_barcode  = 100,
                                                missing_peaks_are_zero = TRUE,
                                                ... ) {
   
   if (sce_save_path == "" | sce_save_path == "." | sce_save_path == "./") {
      stop("Specify some output name to write the sce hdf data, with ",
           "either  'sce_save_path' or 'batch' parameter. ",
           "A directory will be created there.")
   }
   
   # Get list of peaks from the peaks .gtf info file. 
   # Because no all files/cells will have all peaks, need to have a consistant list
   # anything not seen is assumed 0.
   # Obviously needs to be same file for files to merge, so order is the same.
   peak_info <- S4Vectors::DataFrame(read_polyA_peak_file_gtf(peak_info_file) %>%
                                        dplyr::select(peak, chr, start, end, strand, peakdepth, peakgene, misprime)) 
   rownames(peak_info) <- peak_info$peak
   
   
   # Not useing read_tsv because of compression.(segfault.) 
   # NB: https://github.com/tidyverse/readr/issues/610
   the_counts <- read.table(counts_file, sep="\t", header=TRUE, 
                            as.is=TRUE, 
                            colClasses = c("factor","factor", "numeric"),
                            nrows=n_max)
   colnames(the_counts) <- c("peak","cell","count")
   
   # cellprefix to make cell barcodes unique across multi-sample experiments.
   # Batch can be long, but cell_prefix should be short (and unique)
   # Default is batch though.
   # Will also want the barcode (assuming cells named by uniq barcode) later
   cell_to_barcode        <- as.character(the_counts$cell)
   names(cell_to_barcode) <- paste0(cell_prefix, the_counts$cell)
   the_counts$cell        <- factor(names(cell_to_barcode))
   
   
   
   # First construct the sparse matrix of all the counts
   counts_matrix <- Matrix::sparseMatrix(
      i = as.integer(the_counts$peak),
      j = as.integer(the_counts$cell),
      x = as.integer(the_counts$count))
   colnames(counts_matrix) <- levels(the_counts$cell)
   rownames(counts_matrix) <- levels(the_counts$peak)   
   
   # Assume that any gene not present has counts of 0, and match order.
   # (This is used to make consistant gene-list sces to merge later)
   if (missing_peaks_are_zero) {
      
      full_peak_list  <- peak_info$peak 
      absent_genes    <- full_peak_list[ ! full_peak_list %in% rownames(counts_matrix) ]
      stopifnot( all(rownames(counts_matrix) %in% full_peak_list))
      
      counts_matrix.absent <- Matrix::Matrix(0, nrow=length(absent_genes), ncol=ncol(counts_matrix), sparse = TRUE)  
      rownames(counts_matrix.absent) <- absent_genes
      
      counts_matrix <- rbind(counts_matrix, counts_matrix.absent)
      counts_matrix <- counts_matrix[full_peak_list,]
   }
   
   # cell data
   new_colData <- DataFrame(
      cell        = colnames(counts_matrix), 
      barcode     = cell_to_barcode[colnames(counts_matrix)], 
      batch       = batch) 
   
   # Load up the sce.
   sce <- SingleCellExperiment(
      assays = list(counts = counts_matrix ),
      colData = new_colData,
      rowData = DataFrame(peak_info[rownames(counts_matrix),]) )   
   rm(counts_matrix)
   
   # Subset cells with enough reads (since most are rubbish barcodes.)
   sce$total_reads <- DelayedMatrixStats::colSums2(counts(sce))
   
   
   # Save it
   sce <- saveHDF5SummarizedExperiment(sce[,sce$total_reads >= min_reads_per_barcode], 
                                       sce_save_path, ... ) 
  
   return(sce)
}









#' lalalla
#' @export
load_peaks_counts_files_into_sce <- function(
         counts_file_list, batch_names, # Order must match!
         peak_info_file, 
         sce_save_path,
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
   # interium sces need to be saved
   sce_save_paths <- paste0(tempfile(),"_",batch_names)
   sce <- do.call(cbind,
       mapply(FUN=load_peaks_counts_into_sce, 
                       counts_file   = counts_file_list,
                       batch         = batch_names,   
                       sce_save_path = sce_save_paths, 
                       cell_prefix   = batch_names,
                       MoreArgs = list(
                          peak_info_file=peak_info_file,
                          missing_peaks_are_zero = TRUE,
                          #replace is not passed through, should be in tmp
                          ...
      )))

   # Save it (and use replace here.)
   sce <- saveHDF5SummarizedExperiment(sce, sce_save_path, replace=replace)
   
   # Delete the temp files
   file.remove(file.path(sce_save_paths,"assays.h5"))
   file.remove(file.path(sce_save_paths,"se.rds"))
   unlink(sce_save_paths)
   
   return(sce)
}




#' @export
load_peaks_counts_dir_into_sce <- function(
   counts_file_dir, 
   peak_info_file,
   sce_save_path,
   batch_regex = '(.*).tab.gz',
   ... ) {
   
   counts_file_list <- list.files(counts_file_dir)
   batch_names     <-  sub( batch_regex,'\\1', basename(counts_file_list))
   
   sce <- load_peaks_counts_files_into_sce(
      counts_file_list = file.path(counts_file_dir,counts_file_list), 
      batch_names      = batch_names, # Order must match!
      peak_info_file   = peak_info_file, 
      sce_save_path    = sce_save_path,
      ... )
   
   return(sce)
}



   
   





load_peaks_counts_into_sce.viah5 <- function(counts_file, 
                                             peak_info_file, 
                                             batch="",
                                             sce_save_path = batch,
                                             cell_prefix   = batch,
                                             n_max=-1,
                                             min_reads_per_barcode  = 100,
                                             missing_peaks_are_zero = TRUE,
                                             ... ) {
   
   if (sce_save_path == "" | sce_save_path == "." | sce_save_path == "./") {
      stop("Specify some output name to write the sce hdf data, with ",
           "either  'sce_save_path' or 'batch' parameter. ",
           "A directory will be created there.")
   }
   
   # Get list of peaks from the peaks .gtf info file. 
   # Because no all files/cells will have all peaks, need to have a consistant list
   # anything not seen is assumed 0.
   # Obviously needs to be same file for files to merge, so order is the same.
   peak_info <- S4Vectors::DataFrame(read_polyA_peak_file_gtf(peak_info_file) %>%
                                        dplyr::select(peak, chr, start, end, strand, peakdepth, peakgene, misprime)) 
   rownames(peak_info) <- peak_info$peak
   
   
   # Not useing read_tsv because of compression.(segfault.) 
   # NB: https://github.com/tidyverse/readr/issues/610
   the_counts <- read.table(counts_file, sep="\t", header=TRUE, 
                            as.is=TRUE, 
                            colClasses = c("factor","factor", "numeric"),
                            nrows=n_max)
   colnames(the_counts) <- c("peak","cell","count")
   
   # cellprefix to make cell barcodes unique across multi-sample experiments.
   # Batch can be long, but cell_prefix should be short (and unique)
   # Default is batch though.
   # Will also want the barcode (assuming cells named by uniq barcode) later
   cell_to_barcode        <- as.character(the_counts$cell)
   names(cell_to_barcode) <- paste0(cell_prefix, the_counts$cell)
   the_counts$cell        <- factor(names(cell_to_barcode))
   
   
   
   # First construct the sparse matrix of all the counts
   counts_matrix <- Matrix::sparseMatrix(
      i = as.integer(the_counts$peak),
      j = as.integer(the_counts$cell),
      x = as.integer(the_counts$count))
   colnames(counts_matrix) <- levels(the_counts$cell)
   rownames(counts_matrix) <- levels(the_counts$peak)   
   
   # Assume that any gene not present has counts of 0, and match order.
   # (This is used to make consistant gene-list sces to merge later)
   if (missing_peaks_are_zero) {
      
      full_peak_list  <- peak_info$peak 
      absent_genes    <- full_peak_list[ ! full_peak_list %in% rownames(counts_matrix) ]
      stopifnot( all(rownames(counts_matrix) %in% full_peak_list))
      
      counts_matrix.absent <- Matrix::Matrix(0, nrow=length(absent_genes), ncol=ncol(counts_matrix), sparse = TRUE)  
      rownames(counts_matrix.absent) <- absent_genes
      
      counts_matrix <- rbind(counts_matrix, counts_matrix.absent)
      counts_matrix <- counts_matrix[full_peak_list,]
   }
   
   
   
   
   # First save the matrix to hdf5 on disk, filter, ~then~ make the sce
   # Avoids rowsum calculations on full sparse matrix.
   the_temp_h5 <- tempfile()
   counts_matrix.ondisk <- writeTENxMatrix(counts_matrix, the_temp_h5) # rownames truncated??
   
   
   # Up to this point the whole matrix is in memory in sparse format - 
   # best to load a sample at a time with large datasets, then merge hdf5-ly
   # (maybe saveHDF5Array can do this directly, but need to calculate rowsums
   #  regardless - may as well flatten it now.)
   rm(counts_matrix)
   
   if (any(is.na(cell_to_barcode[colnames(counts_matrix.ondisk)]))) {
      stop( paste0(
         "Cell names constructed from cell barcode and ",
         "cell prefix '",cell_prefix, "' do not match sce object cell names. ",
         "This can happen if they are longer than 26 (e.g 16+11) characters combined, ",
         "(saw ", max(nchar(cell_to_barcode))," + ",nchar(cell_prefix),"). ",
         "Consider specifing a smaller character unique cell_prefix/batch (e.g.<11)"))
   }
   
   
   
   # Subset of cells with enough reads (since most are rubbish barcodes.)
   total_readcounts <- DelayedMatrixStats::colSums2(counts_matrix.ondisk)
   new_colData <- DataFrame(
      cell        = colnames(counts_matrix.ondisk), 
      barcode     = cell_to_barcode[colnames(counts_matrix.ondisk)], 
      batch       = batch,
      total_reads = total_readcounts )
   keepers        <- new_colData$total_reads >= min_reads_per_barcode
   
   
   # Load up the sce.
   sce <- SingleCellExperiment(
      assays = list(counts = counts_matrix.ondisk[,keepers] ),
      colData = new_colData[keepers,],
      rowData = DataFrame(peak_info[rownames(counts_matrix.ondisk),]) )
   
   # Save it
   sce <- saveHDF5SummarizedExperiment(sce, sce_save_path, ... ) 
   file.remove( the_temp_h5)
   
   # Sanity check batch+barcode = cell_name
   # Due to error on long batch names/barcode.
   # Longer than 11 + 16 = 26 might cause truncation of the 
   # end of the cell barcode in cell name. Might be version dependant though?
   if (!all(paste0(cell_prefix, sce$barcode) == colnames(sce))) {
      stop( paste0(
         "Cell names constructed from cell barcode and ",
         "cell prefix '",cell_prefix, "' do not match sce object cell names. ",
         "This can happen if they are longer than 26 characters combined, ",
         "(saw ", max(nchar(cell_to_barcode))," + ",nchar(cell_prefix),"). ",
         "Consider specifing a smaller (e.g. <11 character unique cellprefix (or batch)?"))
   }
   
   
   return(sce)
}
