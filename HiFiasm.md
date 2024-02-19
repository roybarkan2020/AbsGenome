# Primery assembly (and QC) with HiFi and Hi-C reads

This is a general workflow for genome assembly using HiFi and Hi-C reads. Start with adapters removal and reads QC. At the end of this workflow, you should have your primary assembly (contigs level, polished with Hi-C reads) and basic stats. With Hifiasm, you can use different modes according to your data availability (i.e., HiFi reads only, Hi-C integration, Trio binning, Ultra-long ONT integration). 

## Convert .bam to .fastq and remove reads with remnant PacBio adapter sequences (https://github.com/sheinasim/HiFiAdapterFilt)
```
hifiadapterfilt.sh -p HiFi_reads.bam -t 16 -o HiFi_reads_filtered.fq.gz
```

## Assembly with HiFi (adapter-free) and Hi-C reads (https://hifiasm.readthedocs.io/en/latest/hic-assembly.html#hic-assembly](https://github.com/chhylp123/hifiasm)
```
hifiasm -o output_name.asm -t16 -h1 HiC_read1.fq.gz -h2 HiC_read2.fq.gz HiFi_reads_filtered.fq.gz
```

## gfa to fasta (hifiasm primery assembly) 
```
awk '/^S/{print ">"$2;print $3}' output_name.hic.p_ctg.gfa > output_name.hic.p_ctg.fa
```

## Run BUSCO on the assembly (https://busco.ezlab.org/busco_userguide.html#getting-started)
```
busco -i output_name.hic.p_ctg.fa -m genome --auto-lineage-euk
```

## Run QUAST to evaluates genome assemblies (https://quast.sourceforge.net/docs/manual.html#sec2)
```
 python quast.py --threads 16 --eukaryote --large -o quast_output_name output_name.hic.p_ctg.fa
```

## Run Meryl ,genomic k-mer counter and sequence utility (https://github.com/marbl/merqury)
```
meryl print \
  union \
    count k=20 HiFi_reads_filtered.fq.gz output.meryldb
```


## Use k-mer counts from Meryl as input for Merqury (https://github.com/marbl/merqury) and GenomeScope https://github.com/schatzlab/genomescope to assess quality and completeness of genome assemblies 

```
merqury.sh output.meryldb output_name.hic.p_ctg.fa output_merqury
```
GenomeScope can be used as command line or in the webtool http://genomescope.org/
```
genomescope.R output.meryldb k-mer_length read_length output_dir
```
