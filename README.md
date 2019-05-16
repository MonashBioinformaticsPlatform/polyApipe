# polyApipe



## polyApipe.py

Count reads in polyA regions with polyApipe.py script on a directory or bam file:

```
polyApipe.py -i demo/ -o demo
polyApipe.py -i demo/SRR5259422_demo.bam -o SRR5259422_demo
```

## polyApiper R package


```
library (polyApiper)

sce <- load_peaks_counts_dir_into_sce("demo_counts/", "demo_polyA_peaks.gff", "./demo_sce") 

...

```





