
document :
	cd polyApiper ; Rscript -e "devtools::document()"

check : document
	cd polyApiper ; rm -f polyApiper_*.tar.gz
	cd polyApiper ; R CMD build .
	cd polyApiper ; R CMD check polyApiper_*.tar.gz
	cd polyApiper ; rm polyApiper_*.tar.gz

site : document
	Rscript -e 'rmarkdown::render("README.md",output_file="index.html",output_dir="docs")'
	Rscript -e 'rmarkdown::render("documentation/polyApipe.Rmd",output_dir="docs")'
	cd polyApiper ; Rscript -e 'pkgdown::build_site(new_process=FALSE)'

