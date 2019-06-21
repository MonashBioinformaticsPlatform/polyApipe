# polyApipe

polyApipe is a pipeline for examining Alternative PolyAdenylation (APA) from 10X Genomics single-cell RNA Sequencing. It consists of a Python 3 part, `polyApipe.py` which performs peak calling and UMI counting, and and R part, `polyApiper`, which analyses the resultant UMI counts.


## Install

Install the Python 3 part, from the command line:

```
pip install git+https://github.com/swbioinf/polyApipe
```

Install the R part in R:

```
BiocManager::install("pfh/weitrix")
BiocManager::install("swbioinf/polyApipe/polyApiper")
```


## polyApipe.py


Count reads in polyA regions with polyApipe.py script on a directory or bam file:

```
polyApipe.py -i demo/ -o demo
polyApipe.py -i demo/SRR5259422_demo.bam -o SRR5259422_demo
```


## polyApiper R package

```
library(polyApiper)

# - Get appropriate ENSEMBL annotations and DNA sequence from AnnotationHub
# - Classify genomic regions into exon,intron,3'UTR,extension
do_ensembl_organism(
    out_path="mouse_ens94", 
    species="Mus musculus", 
    version="94")

# - Load output from polyApipe.py
# - Produce HDF5Array SingleCellExperiment objects containing counts
# - Perform further analysis steps
do_pipeline(
    out_path="demo_banquet", 
    counts_files="demo_counts", 
    peak_info_file="demo_polyA_peaks.gff", 
    organism="mouse_ens94")

# - Load results (individual objects are lazy-loaded when accessed)
organism <- load_banquet("mouse_ens94")
banq <- load_banquet("demo_banquet")
names(banq)
```

### Lower level processing in R

```
sce <- load_peaks_counts_dir_into_sce("demo_counts/", "demo_polyA_peaks.gff", "./demo_sce") 
...

```





