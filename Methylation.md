# Methylation calling from PacBio's HiFi reads 
**This workflow will generate site methylation probabilities from mapped HiFi reads.** 
For additional info about the tool go to [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools)
## Map HiFi reads (raw with 5mC base modification tags and values - MM/ML) to the final genome assembly using minimap2 or pbmm2 (A minimap2 SMRT wrapper for PacBio data - https://github.com/PacificBiosciences/pbmm2)
This step must be done using soft-clipping (-Y flag) to allow the MM and ML tags to remain valid and interpretable. Please note that this step takes time and will produce a large BAM file (proportional to the genome size). 
```
minimap2 -ax map-pb -Y scaffolds.fa HiFi_reads_filtered.fq > CpG.sam
```

## Converting the sam alignment file to a sorted, indexed bam file using samtools
```
samtools view -Sb -o CpG.bam CpG.sam
```
```
samtools sort -O bam -o sorted_CpG.bam CpG.bam
```

## Make sure that the MM/ML (methylation tags) were retained in your output file from the last stap (CpG.bam file)
```
samtools view aln_CpG.bam | head -n 1 | tr '\t' '\n' | less -S
```

## Index BAM
```
samtools index sorted_CpG.bam
```

## pb-CpG-tools
Generate site methylation probabilities from mapped HiFi reads. This will produce two output files - bed and bigwig file. for additional options and output description, go to https://github.com/PacificBiosciences/pb-CpG-tools or type ```aligned_bam_to_cpg_scores --help```
```
aligned_bam_to_cpg_scores \
  --bam sorted_CpG.bam \
  --output-prefix CpG \
  --model models/pileup_calling_model.v1.tflite \
  --threads 16
```

You can now visualise the data with IGV by uploading the following file: 
1. Genome file and index
2. The bam file of the soft-clipped HiFi reads to the genome (the sorted and indexed output of minimap2)
3. The outputs from the aligned_bam_to_cpg_scores step (bed and bigwig files).
Once uploaded, right-click the bam file track and choose "Color alignments by -> base modification (5mC)". 
