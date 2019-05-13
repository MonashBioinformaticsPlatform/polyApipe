This demo dataset is a a biased subset of reads from  SRR5259422 and SRR5259354, 
part of the 10X 1 Million neurons dataset: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.3.0/1M_neurons 
It should take seconds to run.

#commands to generate (enriching for possible polyA read just to keep input sizes down):
samtools view -h SRR5259422_E18_20161004_Neurons_Sample_12.bam | grep -v 114M | grep -E '(^@|AAAAAA|TTTTTT)' | head -n 13034  | samtools view -S -b >  demo/SRR5259422_demo.bam
samtools view -h SRR5259354_E18_20160930_Neurons_Sample_21.bam | grep -v 114M | grep -E '(^@|AAAAAA|TTTTTT)' | head -n 13034  | samtools view -S -b >  demo/SRR5259354_demo.bam
samtools index demo/SRR5259422_demo.bam
samtools index demo/SRR5259354_demo.bam



