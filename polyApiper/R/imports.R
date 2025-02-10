
#'
#' @import methods
#'
#' @import weitrix
#' @import topconfects
#'
#' @importFrom stats cor median runif
#' @importFrom utils data
#'
#' @importFrom rmarkdown render
#'
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join inner_join distinct tibble row_number
#' @importFrom readr write_csv
#' @importFrom tidyr unnest
#' @importFrom purrr map map_lgl map_chr pmap_df
#' @importFrom jsonlite read_json write_json
#' @importFrom assertthat assert_that
#'
#' @importFrom S4Vectors
#'     metadata
#'
#' @importFrom BiocParallel bpparam
#'
#' @importFrom HDF5Array
#'     loadHDF5SummarizedExperiment saveHDF5SummarizedExperiment
#'     quickResaveHDF5SummarizedExperiment
#'     setHDF5DumpDir getHDF5DumpDir
#'
#' @importFrom DelayedArray
#'     getRealizationBackend setRealizationBackend rowsum realize
#'     sweep
#'     rowSums colSums rowMeans colMeans
### Seem to need to import rowSums, colSums, rowMeans, colMeans even though DelayedArray is monkeypatching base so that they are generics or something.
#'
#' @importFrom SingleCellExperiment 
#'     SingleCellExperiment counts logcounts sizeFactors sizeFactors<- reducedDim
#'
#' @importFrom SummarizedExperiment
#'     SummarizedExperiment assay assay<- rowData rowData<- rowRanges colData colData<-
#'
#' @importFrom scran
#'     quickCluster computeSumFactors
#'
#' @importFrom scater
#'     logNormCounts
#'
#' @importFrom scuttle
#'     uniquifyFeatureNames
#'
#' @importFrom BiocGenerics
#'     rbind cbind
## No longer provided by BiocGenerics: rowSums colSums rowMeans colMeans
#'
#' @importFrom plyranges
#'     filter select mutate group_by ungroup arrange summarize
#'     anchor_5p anchor_3p
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
#' @importFrom limma lmFit vooma wilcoxGST
#'
#' @importFrom msigdbr msigdbr
#'
#' @importFrom GSEABase description geneIds
#'
#' @importFrom shiny 
#'      shinyApp fluidPage titlePanel tabPanel navlistPanel
#'      selectizeInput updateSelectizeInput
#'      observeEvent
#'      h1 h2 h3
#'
#' @importFrom DT
#'      DTOutput renderDT dataTableProxy selectRows
#'
NULL


