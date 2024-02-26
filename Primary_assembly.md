# Primary assembly (and QC) with HiFi and Hi-C reads

This is a general workflow for genome assembly of a non-model organism (tropical abalone, Haliotis asinina, marine invertebrate) using PacBio HiFi reads, and Hi-C (Omni-C approach, Dovetail genomics) reads. 
Start with adapters removal and reads QC. At the end of this workflow, you should have your primary assembly (contigs level, polished with Hi-C reads) and basic stats. With Hifiasm, you can use different modes according to your data availability (i.e., HiFi reads only, Hi-C integration, Trio binning, Ultra-long ONT integration). 

## Convert .bam to .fastq and remove reads with remnant PacBio adapter sequences (https://github.com/sheinasim/HiFiAdapterFilt)
Usually, you will get your HiFi reads as bam files. this command will convert the file to fastq format (the required input for Hifiasm) and remove adapter sequences and remnant reads.  
```
hifiadapterfilt.sh -p HiFi_reads.bam -t 16 -o HiFi_reads_filtered.fq.gz
```

## Assembly with HiFi (adapter-free) and Hi-C reads (https://hifiasm.readthedocs.io/en/latest/hic-assembly.html#hic-assembly](https://github.com/chhylp123/hifiasm)
If you have Hi-C (Illumina) paired reads from multiple lanes, you should concatenate/merge all forward (-h1) reads together and reverse (-h2) reads together. 

```
hifiasm -o output_name.asm -t16 -h1 HiC_read1.fq.gz -h2 HiC_read2.fq.gz HiFi_reads_filtered.fq.gz
```

## gfa to fasta (hifiasm primary assembly) 
For the next stages, you will want to convert the assembly files from .gfa files to .fa files. 
```
awk '/^S/{print ">"$2;print $3}' output_name.hic.p_ctg.gfa > output_name.hic.p_ctg.fa
```

## Run BUSCO on the assembly (https://busco.ezlab.org/busco_userguide.html#getting-started)
Run BUSCO with the appropriate lineage database. 
```
busco -i output_name.hic.p_ctg.fa -m genome --auto-lineage-euk
```

## Run QUAST to evaluate genome assemblies (https://quast.sourceforge.net/docs/manual.html#sec2)
```
 python quast.py --threads 16 --eukaryote --large -o quast_output_name output_name.hic.p_ctg.fa
```

## Run Meryl, genomic k-mer counter and sequence utility (https://github.com/marbl/merqury)
```
meryl print \
  union \
    count k=20 HiFi_reads_filtered.fq.gz output.meryldb
```


## Use k-mer counts from Meryl as input for Merqury (https://github.com/marbl/merqury) and GenomeScope/GenomeScope2.0 (https://github.com/schatzlab/genomescope, https://github.com/tbenavi1/genomescope2.0) to assess quality and completeness of genome assemblies 

```
merqury.sh output.meryldb output_name.hic.p_ctg.fa output_merqury
```
GenomeScope can be used as a command line or in the web tool http://genomescope.org/
```
genomescope.R output.meryldb k-mer_length read_length output_dir
```
