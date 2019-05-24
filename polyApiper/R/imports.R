
#'
#' @import methods
#'
#' @import weitrix
#' @import topconfects
#'
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join distinct
#' @importFrom purrr map_lgl
#' @importFrom jsonlite read_json write_json
#' @importFrom assertthat assert_that
#'
#' @importFrom HDF5Array
#'     loadHDF5SummarizedExperiment saveHDF5SummarizedExperiment
#'     quickResaveHDF5SummarizedExperiment
#'     setHDF5DumpDir getHDF5DumpDir
#'
#' @importFrom DelayedArray
#'     getRealizationBackend setRealizationBackend
#'
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment
#'     assay rowData rowData<- rowRanges colData
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
#'     filter select mutate group_by ungroup anchor_5p anchor_3p
#'     join_overlap_inner_directed join_overlap_left_directed
#'     setdiff_ranges_directed
#'     flank_downstream bind_ranges disjoin_ranges_directed
#'     read_gff3 write_gff3
#'
#' @importFrom GenomicRanges trim
#' @importFrom GenomicFeatures genes transcriptsBy threeUTRsByTranscript exonsBy
#'
#' @importFrom GO.db GO.db
#'
#' @importFrom AnnotationHub AnnotationHub
#'
#' @importFrom limma lmFit vooma
#'
NULL


