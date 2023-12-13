# Omni-C analysis (based on code from https://omni-c.readthedocs.io/en/latest/index.html) - Make sure you got the following tools installed: pysam, tabulate, bedtools, matplotlib, pandas, bwa, pairtools, samtools, preseq. 

# Download Omni-C scripts

```{bash}

git clone https://github.com/dovetail-genomics/Omni-C.git

```

# Merge Omni-C reads (fastq.gz):

```{bash}

# Forward reads:
cat AGRF_CAGRF230313856_22FFKJLT3/Abalone_repeat_22FFKJLT3_TGAGCTAG-AACCGTTC_L007_R1.fastq.gz AGRF_CAGRF230313856_22FFKJLT3/Abalone_repeat_22FFKJLT3_TGAGCTAG-AACCGTTC_L008_R1.fastq.gz AGRF_CAGRF230313856_HHNMYDRX3/Abalone_repeat_HHNMYDRX3_TGAGCTAG-GAACGGTT_L001_R1.fastq.gz > omni_merged_all_R1.fastq.gz

# Reverse reads:
cat AGRF_CAGRF230313856_22FFKJLT3/Abalone_repeat_22FFKJLT3_TGAGCTAG-AACCGTTC_L007_R2.fastq.gz AGRF_CAGRF230313856_22FFKJLT3/Abalone_repeat_22FFKJLT3_TGAGCTAG-AACCGTTC_L008_R2.fastq.gz AGRF_CAGRF230313856_HHNMYDRX3/Abalone_repeat_HHNMYDRX3_TGAGCTAG-GAACGGTT_L001_R2.fastq.gz > omni_merged_all_R2.fastq.gz

```
# Index genome file:
```{bash}

singularity run /fast/tmp/containers/samtools-1.16.1.sif samtools faidx Hal_Asi.asm.bp.p_ctg.fa

cut -f1,2 Hal_Asi.asm.bp.p_ctg.fa.fai > hal_asi_ctg.genome

singularity run /fast/tmp/containers/bwa-0.7.17.sif bwa index Hal_Asi.asm.bp.p_ctg.fa
```

# From .fastq to final valid pairs bam file

# Alignment to the genome (long):
```{bash}
singularity run /fast/tmp/containers/bwa-0.7.17.sif bwa mem -5SP -T0 -t16 Genome/Hifiasm/Hal_Asi.asm.bp.p_ctg.fa Genome/OmniC/omni_merged_all_R1.fastq.gz Genome/OmniC/omni_merged_all_R2.fastq.gz -o aligned_omni_all.sam
```
# Recording valid ligation events:
```{bash}
singularity run /fast/tmp/containers/pairtools-1.0.2.sif pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path Genome/Hifiasm/hal_asi_ctg.genome  aligned_omni_all.sam >  parsed.pairsam
```
# Sorting the pairsam file:
```{bash}
singularity run /fast/tmp/containers/pairtools-1.0.2.sif  pairtools sort --nproc 16 --tmpdir=temp/  parsed.pairsam > sorted.pairsam
```
# Removing PCR duplicates
```{bash}
singularity run /fast/tmp/containers/pairtools-1.0.2.sif pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt --output dedup.pairsam sorted.pairsam
```
# Generating .pairs and bam files
```{bash}
singularity run /fast/tmp/containers/pairtools-1.0.2.sif pairtools split --nproc-in 8 --nproc-out 8 --output-pairs mapped.pairs --output-sam unsorted.bam dedup.pairsam
```
# Generating the final bam file (and index)
```{bash}
singularity run /fast/tmp/containers/samtools-1.16.1.sif samtools sort -@16 -T temp/temp.bam -o mapped.PT.bam unsorted.bam

singularity run /fast/tmp/containers/samtools-1.16.1.sif samtools index mapped.PT.bam

```
# Scaffolding using yahs
```{bash}
/home/jc748673/yahs/yahs /home/jc748673/Genome/Hifiasm/Hal_Asi.asm.bp.p_ctg.fa /home/jc748673/Omni-C/Output_omni/mapped.PT.bam
```

# From .pairs to .hic contact matrix
```{bash}
java -Xmx48000m  -Djava.awt.headless=true -jar Omni-C/juicertools.jar pre --threads 16 mapped.pairs contact_map.hic Genome/Hifiasm/hal_asi_ctg.genome
```

# index scaffolds_final.fa
```{bash}
samtools faidx yahs.out_scaffolds_final.fa
```

# create *.chrom.sizes file
```{bash}
awk '{print $1 " " $2}' yahs.out_scaffolds_final.fa.fai > scaffolds_final.chrom.sizes
```

# Generate HiC contact maps
```{bash}
/home/jc748673/Scaffolding/yahs/juicer pre /home/jc748673/Scaffolding/yahs/yahs.out.bin /home/jc748673/Scaffolding/yahs/yahs.out_scaffolds_final.agp /home/jc748673/Genome/Hifiasm/Hal_Asi.asm.bp.p_ctg.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > alignments_sorted.txt.part && mv alignments_sorted.txt.part alignments_sorted.txt
```

# index scaffolds_final.fa
```{bash}
samtools faidx yahs.out_scaffolds_final.fa
```
# create *.chrom.sizes file
```{bash}
awk '{print $1 " " $2}' yahs.out_scaffolds_final.fa.fai > scaffolds_final.chrom.sizes

java -Xmx48000m  -Djava.awt.headless=true -jar /home/jc748673/Scaffolding/Omni-C/juicertools.jar pre /home/jc748673/Scaffolding/yahs/alignments_sorted.txt out.hic.part scaffolds_final.chrom.sizes && mv out.hic.part out.hic
```

# JBAT files for manual curation 
```{bash}
/home/jc748673/Scaffolding/yahs/juicer pre -a -o out_JBAT /home/jc748673/Scaffolding/yahs/yahs.out.bin /home/jc748673/Scaffolding/yahs/yahs.out_scaffolds_final.agp /home/jc748673/Genome/Hifiasm/Hal_Asi.asm.bp.p_ctg.fa.fai >out_JBAT.log 2>&1


(java -jar -Xmx32G /home/jc748673/Scaffolding/Omni-C/juicertools.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)

```


