# polyApipe

polyApipe is a pipeline for examining Alternative PolyAdenylation (APA) from 10X Genomics single-cell RNA Sequencing. It consists of a Python 3 part, `polyApipe.py` which performs peak calling and UMI counting, and and R part, `polyApiper`, which analyses the resultant UMI counts.

* [GitHub repository](https://github.com/MonashBioinformaticsPlatform/polyApipe)
* [polyApipe.py usage example and documentation](https://monashbioinformaticsplatform.github.io/polyApipe/polyApipe.html)
* [R function documentation](https://monashbioinformaticsplatform.github.io/polyApipe/polyApiper/reference/index.html)

<!-- An example analysis is [here](https://bioinformatics.erc.monash.edu/home/sarah.williams/projects/TraudeBeliharz_SingleCell3pEnds/polyApipe_doco.html) -->

[This poster on polyApipe](https://doi.org/10.7490/f1000research.1117076.1) was presented at Oz Single Cell 2019.


## Install dependencies

polyApipe requires various tools:

 * [samtools](http://www.htslib.org/) - [cite](https://www.ncbi.nlm.nih.gov/pubmed/19505943)
 * [UMI-tools](https://github.com/CGATOxford/UMI-tools) - [cite](https://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract)
 * [featureCounts](http://subread.sourceforge.net/) (version 1.5.3 or above) from Subread package - [cite](https://www.ncbi.nlm.nih.gov/pubmed/24227677)
 * [pysam](https://github.com/pysam-developers/pysam) python module

The quickest way to install these is to use [conda](https://docs.conda.io/en/latest/). 
After installing conda, the following will create a suitable environment, 
which can be loaded with `conda activate`. Ensure it is python3.

``` 
conda create -n polyApipe_env  --override-channels -c bioconda -c conda-forge -c anaconda umi_tools=1.0.0-0 pysam samtools subread 
conda activate polyApipe_env 
```

The R part of polyApipe requires Bioconductor. In R:

```
install.packages("BiocManager")
```

## Install

Download polyApipe:

```
git clone https://github.com/MonashBioinformaticsPlatform/polyApipe
```

Test that you can run the Python part:

```
polyApipe/polyApipe.py
```

Install the R part in R. weitrix development happens in step with polyApipe, so install it from git as well.

```
BiocManager::install("pfh/weitrix")
BiocManager::install("MonashBioinformaticsPlatform/polyApipe/polyApiper")
```


## polyApipe.py Python script


Count reads in poly(A) regions with polyApipe.py script on a directory or bam file:

```
polyApipe/polyApipe.py -i demo/ -o demo
polyApipe/polyApipe.py -i demo/SRR5259422_demo.bam -o SRR5259422_demo
```


## polyApiper R package

Process the output of `polyApipe.py` in R.

```
library(polyApiper)


# - Get appropriate ENSEMBL annotations and DNA sequence from AnnotationHub
# - Classify genomic regions into exon,intron,3'UTR,extension
do_ensembl_organism(
    out_path="mouse_ens100", 
    species="Mus musculus", 
    version="100")


# You will need to supply a set of cells to use.
# - Choose these using conventional scRNA-Seq cell filtering methods.
# - Also provide a cell naming function that combines the batch name and barcode
#   consistent with your existing cell naming.

cells_to_use <- c("SRR6129051_AAACCTGAGAGGGCTT", "SRR6129051_AAACCTGCATACGCCG", ...)

cell_name_func <- function(batch,barcode) paste0(batch,"_",barcode)


# Run the R part of the pipeline
# - Load output from polyApipe.py
# - Produce HDF5Array SingleCellExperiment objects containing counts
# - Perform further analysis steps
do_pipeline(
    out_path="demo_banquet", 
    counts_file_dir="demo_counts", 
    peak_info_file="demo_polyA_peaks.gff", 
    organism="mouse_ens100",
    cells_to_use=cells_to_use,
    cell_name_func=cell_name_func)


# Load results (individual objects are lazy-loaded when accessed)
organism <- load_banquet("mouse_ens100")
banq <- load_banquet("demo_banquet")
names(banq)
```


### Finding the most recent available ENSEMBL annotation version

To work out the available ENSEMBL versions, use this R code, substituting in the appropriate species:

```
library(AnnotationHub)
ah <- AnnotationHub()
subset(ah, rdataclass == "EnsDb" & species == "Homo sapiens")$title
subset(ah, rdataclass == "EnsDb" & species == "Mus musculus")$title
```


### Lower level processing in R

Peak counts can be loaded into a SingleCellExperiment object without any further interpretation using:

```
sce <- load_peaks_counts_dir_into_sce("demo_counts/", "demo_polyA_peaks.gff", "./demo_sce") 
```


### BiocParallel

polyApiper uses uses BiocParallel's default parallel processing as its own default, as returned by `BiocParallel::bpparam()`.

If polyApiper hangs or produces strange errors, you can try running it with serial processing:

```
BiocParallel::register( BiocParallel::SerialParam() )
```




