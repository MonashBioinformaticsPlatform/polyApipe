
#'
#' @import methods
#'
#' @import weitrix
#' @import topconfects
#'
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom magrittr %>%
#' @importFrom future %<-%
#' @importFrom dplyr left_join
#' @importFrom purrr map_lgl
#'
#' @importFrom HDF5Array
#'     loadHDF5SummarizedExperiment saveHDF5SummarizedExperiment setHDF5DumpDir getHDF5DumpDir
#'
#' @importFrom DelayedArray
#'     getRealizationBackend setRealizationBackend
#'
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment assay<-
#'
#' @importFrom scran
#'     computeSumFactors
#'
#' @importFrom scater
#'     normalize
#'
#' @importFrom BiocGenerics
#'     rowSums colSums rbind cbind
#'
#' @importFrom plyranges
#'     filter select mutate group_by n anchor_5p anchor_3p join_overlap_inner_directed setdiff_ranges_directed
#'     flank_downstream bind_ranges disjoin_ranges_directed
#'
#' @importFrom GenomicRanges trim
#' @importFrom GenomicFeatures genes transcriptsBy threeUTRsByTranscript exonsBy
#'
NULL
