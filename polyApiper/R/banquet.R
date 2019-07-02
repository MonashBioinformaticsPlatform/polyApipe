
ensure_dir <- function(path) {
    if (!file.exists(path)) {
        ensure_dir(dirname(path))
        dir.create(path)
    }

    stopifnot(file.info(path)$isdir)
}

clear_dump_dir <- function(path) {
    if (!file.exists(path)) return()

    for(file_path in list.files(path, pattern="auto.*\\.h5", full.names=TRUE))
        file.remove(file_path)
    file.remove(path)
}

#' @export
working_in <- function(path, task) {
    ensure_dir(path)
    dump_dir <- file.path(path, "__working__")
    old_dump_dir <- getHDF5DumpDir()
    old_realization <- getRealizationBackend()

    setHDF5DumpDir(dump_dir)
    setRealizationBackend("HDF5Array")
    on.exit({
        setRealizationBackend(old_realization)
        setHDF5DumpDir(old_dump_dir)
        clear_dump_dir(dump_dir)
    })

    task
}

identify <- function(path) {
    bn <- basename(path)

    if (file.info(path)$isdir) {
        what <- "dir"
        if (file.exists(file.path(path,"assays.h5")))
            what <- "hdf5se"
        return(list(path=path, name=bn, what=what))
    }

    ext <- tolower(file_ext(bn))
    name <- file_path_sans_ext(bn)
    if (ext %in% c("gz","bgz")) {
        ext <- tolower(file_ext(name))
        name <- file_path_sans_ext(name)
    }

    what <- "unknown"
    if (ext == "rds") what <- "rds"
    if (ext == "gff3") what <- "gff3"
    if (ext == "txt") what <- "txt"

    list(
        path=path,
        name=name,
        what=what)
}

#' Load directories lazily, and load RDS and HDF5Array summarized experiments
#'
#' Loads a directory as an environment that will recursively call load_banquet
#' on its contents. Also knows how to load .rds files and HDF5Array summarized
#' experiments.
#'
#' If \code{path} is not a character vector, it is returned unchanged, so this
#' function can be used to ensure an object is loaded.
#'
#' @param path A directory containing .rds files, HDF5Array summarized
#'   experiments, and banquets.
#'
#' @export
load_banquet <- function(path) {
    if (!is(path,"character"))
        return(path)

    if (!file.exists(path))
        stop("No such file: ", path)

    ident <- identify(path)

    if (ident$what == "rds")
        return(readRDS(path))

    if (ident$what == "hdf5se")
        return(loadHDF5SummarizedExperiment(path))

    if (ident$what == "gff3")
        return(read_gff3(path))
        #TODO: set genome_info
    
    if (ident$what == "txt")
        return(readLines(path))

    if (ident$what != "dir")
        stop("Don't know how to load: ",path)

    result <- new.env()

    # Load configuration
    config <- list()
    json_path <- file.path(path, "config.json")
    if (file.exists(json_path))
        config <- read_json(json_path)

    # Load specified AnnotationHub resources
    for(name in names(config$AnnotationHub)) {
        ahid <- config$AnnotationHub[[name]]
        eval(substitute(
            delayedAssign(name, get_ah(ahid), assign.env=result),
            list(ahid=ahid)))
    }

    # Load files
    for(file_path in list.files(path, full.names=TRUE)) {
        sub_ident <- identify(file_path)
        if (sub_ident$what != "unknown") {
            eval(substitute(
                delayedAssign(
                    sub_ident$name, load_banquet(file_path), assign.env=result),
                list(file_path=file_path)))
        }
    }

    result
}
