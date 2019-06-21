
document :
	Rscript -e "devtools::document()"

check : document
	R CMD build .
	R CMD check polyApiper_*.tar.gz
	rm polyApiper_*.tar.gz
