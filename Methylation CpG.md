# Methylation calling from PacBio's HiFi reads (https://github.com/PacificBiosciences/pb-CpG-tools)
This workflow will generate site methylation probabilities from mapped HiFi reads. 

## Map HiFi reads (raw with 5mC base modification tags and values - MM/ML) to the genome assembly using minimap2 or pbmm2 (A minimap2 SMRT wrapper for PacBio data)
This step must be done using  soft-clipping (-Y flag) to allow the MM and ML tags to remain valid and interpretable. Please note that this step takes time and will produce a large BAM file (proportional to the genome size). 
```
minimap2 -ax map-pb -Y scaffolds_final.fa HiFi.fastq > aln_CpG.sam
```
## Converting the sam alignment file to a sorted, indexed bam file using samtools
```
samtools view -Sb -o aln_CpG.bam aln_CpG.sam
```
```
samtools sort -O bam -o sorted_aln_CpG.bam aln_CpG.bam
```

## Index BAM
```
samtools index sorted_aln_CpG.bam
```
## pb-CpG-tools
```
pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
  --bam /home/jc748673/sorted_aln_CpG.bam \
  --output-prefix Hal.CpG.pbmm2 \
  --model pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
  --threads 8
```
