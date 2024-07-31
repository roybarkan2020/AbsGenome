Iso-Seq analysis (based on https://ucdavis-bioinformatics-training.github.io/2020-september-isoseq/liz/bioconda/2-bioconda)
Load PacBio (bioconda) tools
module load pbbioconda/20231004 #for HPC

Merge iso-seq flnc from all samples/tissues
ls H_liver-flnc.bam H_gonad-flnc.bam H_eyes-flnc.bam H_gills-flnc.bam H_tentacle-flnc.bam > flnc_all.fofn

isoseq3 cluster flnc_all.fofn clustered.bam --use-qvs

Align Full Length Non Chemiric (FLNC) reads to ref genome
pbmm2 align yahs.out_scaffolds_final.fa clustered.hq.bam aln.all_hq.bam --sort --preset ISOSEQ --log-level INFO

Collapse sorted aligned isoseq (bam)
isoseq3 collapse aln.all.bam alz.collapsed.all.gff

# map align isoseq high-quality isoforms to genome

```{bash}

minimap2 -ax splice -uf -C5  yahs_polished.out_scaffolds_final.fa High_Quality_Isoforms.fasta > aln_iso_hq_polished.sam

```

# sort sam

```{bash}

samtools sort aln_iso_hq_polished.sam -O sam -o aln_iso_hq_polished_sort.sam

```

# https://github.com/GenomeRIK/tama/wiki/TAMA-GO:-ORF-and-NMD-predictions (link to github manual TAMA)

# TAMA collapse not-capped (the BED files can be converted into GFF with exons and mRNA/transcript tags - done here using BED_to_GFF in Galaxy)

```{bash}

python tama_collapse.py -s aln_iso_hq_polished_sort.sam -f yahs_polished.out_scaffolds_final.fa -p iso_hq_tama_polished -x no_cap

```

# TAMA collapse capped (the BED files can be converted into GFF with exons and mRNA/transcript tags - done here using BED_to_GFF in Galaxy)

```{bash}

python tama_collapse.py -s aln_iso_hq_polished_sort.sam -f yahs_polished.out_scaffolds_final.fa -p iso_hq_tama_polished_capped -x capped

```

# TAMA_degradation_signature

```{bash}

python tama_degradation_signature.py -c iso_hq_tama_polished_capped_trans_read.bed -nc iso_hq_tama_polished_trans_read.bed -o tama_degradation_signature

```

# tama merge bed files annotations

```{bash}

python tama_merge.py -f filelist.txt -p merged_annos

```

# filelist.txt for merging beds:

```{bash}

iso_hq_tama_polished_capped.bed	capped	1,1,1	caplib
iso_hq_tama_polished.bed	capped	2,1,1	nocaplib

```

# Converting the bed files into fasta files

```{bash}

bedtools getfasta -name -split -s -fi yahs_polished.out_scaffolds_final.fa -bed isoseq_polished/merged_annos.bed -fo merge_anno_bed_to_fa

```

# Getting open read frames (ORF) from the transcript sequences

```{bash}

python tama_orf_seeker.py -f merge_anno_bed_to_fa -o orf_aa_merged_polished

```

# Split orf_aa_merged_polished into files for blastp (https://github.com/GenomeRIK/tama/wiki/TAMA-GO:-Split-Files). The number of outputs is automatically chosen by tama_fasta_splitter or can be selected by the user. In this case, split into 9 output files. 

```{bash}

python tama_fasta_splitter.py orf_aa_merged_polished split_orf_aa_merged_polished num_splits

```
# Make the Uniprot/Uniref/Whatever protein database file to a blastp database file

```{bash}

makeblastdb -in uniprot_sprot.fasta -dbtype prot

```


# Blasting the amino acid sequences against the Uniprot/Uniref (in this case it was done with UniprotKB_SwissProt) protein database (run for each of the nine  spliced files)

```{bash}

blastp -evalue 1e-10 -num_threads 16 -db uniprot_sprot.fasta -query split_orf_aa_merged_polished.fa

```

# Concatanate all of the nine blastp outputs

```{bash}

cat blastp_1 blastp_2 blastp_3 blastp_4 blastp_5 blastp_6 blastp_7 blastp_8 blastp_9 > blastp_orf

```

# Parsing the Blastp output file for top hits

```{bash}

python tama_orf_blastp_parser.py -b blastp_orf -o blastp_orf_parsing

```

# Create new bed file with CDS regions

python tama_cds_regions_bed_add.py -p ${orf} -a ${bed} -f ${fasta} -o ${outfile}

