# polyApipe



## polyApipe.py

Count reads in polyA regions with polyApipe.py script on a directory or bam file:

```
polyApipe.py -i demo/ -o demo
polyApipe.py -i demo/SRR5259422_demo.bam -o SRR5259422_demo
```

## polyApiper R package

```
# Install from GitHub, with dependencies
BiocManager::install("pfh/weitrix")
BiocManager::install("swbioinfo/polyApipe/polyApiper")

# Load
library(polyApiper)

# Automated pipeline

# - Get appropriate ENSEMBL annotations and DNA sequence from AnnotationHub
# - Classify genomic regions into exon,intron,3'UTR,extension
do_ensembl_organism("mouse_ens94", "Mus musculus", "94")

# - Load output from polyApipe.py
# - Produce HDF5Array SingleCellExperiment objects containing counts
# - Perform further analysis steps
do_pipeline("demo_banquet", "demo_counts", "demo_polyA_peaks.gff", "mouse_ens94")

# - Load results (individual objects are lazy-loaded when accessed)
organism <- load_banquet("mouse_ens94")
banq <- load_banquet("demo_banquet")
names(banq)


# Lower level processing
sce <- load_peaks_counts_dir_into_sce("demo_counts/", "demo_polyA_peaks.gff", "./demo_sce") 

...

```





