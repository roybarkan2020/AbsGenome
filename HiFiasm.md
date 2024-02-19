# Convert .bam to .fastq and remove reads with remnant PacBio adapter sequences (https://github.com/sheinasim/HiFiAdapterFilt)
```
hifiadapterfilt.sh -p HiFi_reads.bam -t 16 -o HiFi_reads_filtered.fq.gz
```

# Assembly with HiFi (adapter-free) and Hi-C reads (https://hifiasm.readthedocs.io/en/latest/hic-assembly.html#hic-assembly](https://github.com/chhylp123/hifiasm)
```
hifiasm -o output_name.asm -t16 -h1 HiC_read1.fq.gz -h2 HiC_read2.fq.gz HiFi_reads_filtered.fq.gz
```

# gfa to fasta (hifiasm primery assembly) 
```
awk '/^S/{print ">"$2;print $3}' output_name.hic.p_ctg.gfa > output_name.hic.p_ctg.fa
```

# Run BUSCO on the assembly (https://busco.ezlab.org/busco_userguide.html#getting-started)
```
busco -i output_name.hic.p_ctg.fa -m genome --auto-lineage-euk
```

# Run QUAST to evaluates genome assemblies (https://quast.sourceforge.net/docs/manual.html#sec2)
```
 python quast.py --threads 16 --eukaryote --large -o quast_output_name output_name.hic.p_ctg.fa
```

# Run Meryl ,genomic k-mer counter and sequence utility, and Merqury ,reference-free quality, completeness, and phasing assessment for genome assemblies (https://github.com/marbl/merqury)
```
meryl print \
  union \
    count k=20 HiFi_reads_filtered.fq.gz output.meryldb
```
```
merqury.sh output.meryldb output_name.hic.p_ctg.fa output_merqury
