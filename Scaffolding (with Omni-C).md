# Hi-C analysis 

The following provides the steps used for Scaffolding the abalone genome using Omni-C data. 
**Omni-C approach**, code for the analysis and Hi-C maps is based on workflow from [Omni-C’s documentation - Dovetail Genomics](https://omni-c.readthedocs.io/en/latest/index.html). The scaffolding was done using [YaHS](https://github.com/c-zhou/yahs)

## Download Omni-C scripts

```
git clone https://github.com/dovetail-genomics/Omni-C.git
```

## Merge Omni-C reads (fastq.gz):
In some cases, your paired-end Hi-C data will be sequenced in multiple lanes. Usually, you will receive the data for each lane separately. In this case, you should merge all forward reads into one file and all reverse reads into one file.  To do so, use:

1. Forward reads:
```
cat R1_lane001.fastq.gz R1_lane002.fastq.gz R1_lane003.fastq.gz > OmniC.R1.fastq.gz
```
2. Reverse reads:
```
cat R2_lane001.fastq.gz R2_lane002.fastq.gz R2_lane003.fastq.gz > OmniC.R2.fastq.gz
```

## Pre-alignment steps 

### Index genome assembly file:

```
samtools faidx genome.fa
```

### Use the index file to generate the genome file

```
cut -f1,2 genome.fa.fai > genome.genome
```

### Generate a BWA index file

```
bwa index genome.fa
```

Library QC and library complexity were both performed by the sequencing service provider as part of the project. **These steps are not described in this workflow.** 
For the **library QC** step [click here](https://omni-c.readthedocs.io/en/latest/library_qc.html#)
For the **library complexity** step [click here](https://omni-c.readthedocs.io/en/latest/library_qc.html#library-complexity) 

## From .fastq to final valid pairs bam file

### Alignment to the genome (long):
```
bwa mem -5SP -T0 -t16 genome.fa OmniC.R1.fastq.gz OmniC.R2.fastq.gz -o aligned.sam
```
### Recording valid ligation events:
```
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path genome.genome  aligned.sam >  parsed.pairsam
```
### Sorting the pairsam file:
```
pairtools sort --nproc 16 --tmpdir=temp/  parsed.pairsam > sorted.pairsam
```
### Removing PCR duplicates
```
pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt --output dedup.pairsam sorted.pairsam
```
### Generating .pairs and bam files
```
pairtools split --nproc-in 8 --nproc-out 8 --output-pairs mapped.pairs --output-sam unsorted.bam dedup.pairsam
```
### Generating the final bam file (and index)

```
samtools sort -@16 -T temp/temp.bam -o mapped.bam unsorted.bam
```

```
samtools index mapped.bam
```

<table><tr><td>There is an option of running it all using one command by piping all of the above. This can be done using: </td></tr></table>

```
bwa mem -5SP -T0 -t<cores> <ref.fa> <OmniC.R1.fastq.gz> <OmniC.R2.fastq.gz>| \
pairtools parse --min-mapq 40 --walks-policy 5unique \
--max-inter-align-gap 30 --nproc-in <cores> --nproc-out <cores> --chroms-path <ref.genome> | \
pairtools sort --tmpdir=<full_path/to/tmpdir> --nproc <cores>|pairtools dedup --nproc-in <cores> \
--nproc-out <cores> --mark-dups --output-stats <stats.txt>|pairtools split --nproc-in <cores> \
--nproc-out <cores> --output-pairs <mapped.pairs> --output-sam -|samtools view -bS -@<cores> | \
samtools sort -@<cores> -o <mapped.PT.bam>;samtools index <mapped.PT.bam>
```

<ins>**Moving on to the scaffolding:** </ins>

## Scaffolding using [YaHS](https://github.com/c-zhou/yahs)

```
yahs genome.fa mapped.bam -o yahs.out
```

### From .pairs to .hic contact matrix


This step is done using [Juicer - AidenLab](https://github.com/aidenlab/juicer). 
More tools from **Aiden Lab** are available for manual curation and analysis of Hi-C data [AidenLab](https://github.com/aidenlab)

```
java -Xmx48000m  -Djava.awt.headless=true -jar juicertools.jar pre --threads 16 mapped.pairs contact_map.hic genome.genome
```

### index scaffolds_final.fa

```
samtools faidx yahs.out.fa
```

### create *.chrom.sizes file

```
awk '{print $1 " " $2}' yahs.out.fa.fai > yahs_scaffolds_final.chrom.sizes
```

### Generate HiC contact maps

```
juicer pre yahs.out.bin yahs.out.agp genome.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > alignments_sorted.txt.part && mv alignments_sorted.txt.part alignments_sorted.txt
```

```
java -Xmx48000m  -Djava.awt.headless=true -jar juicertools.jar pre alignments_sorted.txt out.hic.part yahs_scaffolds_final.chrom.sizes && mv out.hic.part out.hic
```

out.hic can be used to visualise the Hi-C heatmap using **JuiceBox**



