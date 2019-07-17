

#' @export
do_pipeline_report <- function(
        out_path, 
        input=out_path, 
        title="polyApiper results overview") {

    input <- load_banquet(input)
    wd <- getwd()
    
    ensure_dir(out_path)
    
    rmarkdown::render(
        system.file("rmd","pipeline_overview.Rmd", package="polyApiper"),
        output_dir=out_path,
        output_file=file.path(out_path,"overview.html"),
        params=list(
            input=input,
            title=title,
            wd=wd))
    
    invisible(NULL)
}