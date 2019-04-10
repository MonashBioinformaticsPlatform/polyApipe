
#' Turn a directory into an environment of lazy-loading objects
#'
#' @param path A directory containing .rds files, HDF5Array summarized experiments, and banquets.
#'
#' @export
load_banquet <- function(path) {
    if (grepl("\\.rds$", path, ignore.case=TRUE))
        return(readRDS(path))

    if (!file.info(path)$isdir)
        stop("Don't know how to load: ",path)

    if (file.exists(file.path(path,"assays.h5")))
        return(loadHDF5SummarizedExperiment(path))

    result <- new.env()

    for(filename in list.files(path)) {
        name <- sub("\\.rds$","", filename)
        file_path <- file.path(path, filename)
        result[[name]] %<-% load_banquet(file_path)
    }

    result
}
