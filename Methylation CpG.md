# Map HiFi reads (raw with CpG flags) to the genome assembly using minimap2 (soft clipping mode -Y)

minimap2 -ax map-pb -Y /home/jc748673/Scaffolding/yahs/yahs.out_scaffolds_final.fa /home/jc748673/Genome/Hal_Asi_HiFi.fastq > aln_CpG.sam

# Converting the sam alignment file to a sorted, indexed bam file using samtools

samtools view -Sb -o aln_CpG.bam aln_CpG.sam

samtools sort -O bam -o sorted_aln_CpG.bam aln_CpG.bam

# Index BAM

samtools index sorted_aln_CpG.bam

# pb-CpG-tools

pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
  --bam /home/jc748673/sorted_aln_CpG.bam \
  --output-prefix Hal.CpG.pbmm2 \
  --model pb-CpG-tools-v2.3.1-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
  --threads 8
