# polyApipe

polyApipe is a pipeline for examining Alternative PolyAdenylation (APA) from 10X Genomics single-cell RNA Sequencing. It consists of a Python 3 part, `polyApipe.py` which performs peak calling and UMI counting, and and R part, `polyApiper`, which analyses the resultant UMI counts.

An example analysis is [here](https://bioinformatics.erc.monash.edu/home/sarah.williams/projects/TraudeBeliharz_SingleCell3pEnds/polyApipe_doco.html)

[This poster on polyApipe](https://doi.org/10.7490/f1000research.1117076.1) was presented at Oz Single Cell 2019.


## Install

Install the Python 3 part, from the command line:

```
pip install --user git+https://github.com/MonashBioinformaticsPlatform/polyApipe
```

You many need to add the polyApipe.py script to your \$PATH e.g. `export PATH=$PATH:/home/myusername/.local/bin/`

Then, install the required tools:

 * [samtools](http://www.htslib.org/) - [cite](https://www.ncbi.nlm.nih.gov/pubmed/19505943)
 * [UMI-tools](https://github.com/CGATOxford/UMI-tools) - [cite](https://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract)
 * [featureCounts](http://subread.sourceforge.net/) (version 1.5.3 or above) from Subread package - [cite](https://www.ncbi.nlm.nih.gov/pubmed/24227677)
 * [pysam](https://github.com/pysam-developers/pysam) python module

The quickest way is to use [conda](https://docs.conda.io/en/latest/). 
After installing conda, the following will create a suitable environment, 
which can be loaded with `conda activate`. Ensure it is python3.

``` 
conda create -n polyApipe_env  --override-channels -c bioconda -c conda-forge -c anaconda umi_tools=1.0.0-0 pysam samtools subread 
conda activate polyApipe_env 
```

Install the R part in R:

```
BiocManager::install("pfh/weitrix")
BiocManager::install("MonashBioinformaticsPlatform/polyApipe/polyApiper")
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






