# map align isoseq high-quality isoforms to genome

```{bash}

minimap2 -ax splice -uf -C5  /home/jc748673/Scaffolding/test_polish/yahs_polished.out_scaffolds_final.fa /home/jc748673/Genome/IsoSeq/IsoSeq_Analysis/High_Quality_Isoforms.fasta > aln_iso_hq_polished.sam

```

# sort sam

```{bash}

samtools sort aln_iso_hq_polished.sam -O sam -o aln_iso_hq_polished_sort.sam

```

# https://github.com/GenomeRIK/tama/wiki/TAMA-GO:-ORF-and-NMD-predictions (link to github manual TAMA)

# TAMA collapse not-capped (the BED files can be converted into GFF with exons and mRNA/transcript tags - done here using BED_to_GFF in Galaxy)

```{bash}

singularity run /fast/tmp/containers/python-2.7.17.sif python /home/jc748673/tama/tama_collapse.py -s /home/jc748673/Scaffolding/test_polish/isoseq_polished/aln_iso_hq_polished_sort.sam -f /home/jc748673/Scaffolding/test_polish/yahs_polished.out_scaffolds_final.fa -p /home/jc748673/Scaffolding/test_polish/isoseq_polished/iso_hq_tama_polished -x no_cap

```

# TAMA collapse capped (the BED files can be converted into GFF with exons and mRNA/transcript tags - done here using BED_to_GFF in Galaxy)

```{bash}

singularity run /fast/tmp/containers/python-2.7.17.sif python /home/jc748673/tama/tama_collapse.py -s /home/jc748673/Scaffolding/test_polish/isoseq_polished/aln_iso_hq_polished_sort.sam -f /home/jc748673/Scaffolding/test_polish/yahs_polished.out_scaffolds_final.fa -p /home/jc748673/Scaffolding/test_polish/isoseq_polished/iso_hq_tama_polished_capped -x capped

```

# TAMA_degradation_signature

```{bash}

singularity run /fast/tmp/containers/python-2.7.17.sif python /home/jc748673/tama/tama_go/file_stats/tama_degradation_signature.py -c /home/jc748673/Scaffolding/test_polish/isoseq_polished/iso_hq_tama_polished_capped_trans_read.bed -nc /home/jc748673/Scaffolding/test_polish/isoseq_polished/iso_hq_tama_polished_trans_read.bed -o /home/jc748673/Scaffolding/test_polish/isoseq_polished/tama_degradation_signature

```

# tama merge bed files annotations

```{bash}

python /home/jc748673/tama/tama_merge.py -f /home/jc748673/Scaffolding/test_polish/isoseq_polished/filelist.txt -p /home/jc748673/Scaffolding/test_polish/isoseq_polished/merged_annos

```

# filelist.txt for merging beds:

```{bash}

/home/jc748673/Scaffolding/test_polish/isoseq_polished/iso_hq_tama_polished_capped.bed	capped	1,1,1	caplib
/home/jc748673/Scaffolding/test_polish/isoseq_polished/iso_hq_tama_polished.bed	capped	2,1,1	nocaplib

```

# Converting the bed files into fasta files

```{bash}

bedtools getfasta -name -split -s -fi /home/jc748673/Scaffolding/test_polish/yahs_polished.out_scaffolds_final.fa -bed /home/jc748673/Scaffolding/test_polish/isoseq_polished/merged_annos.bed -fo /home/jc748673/Scaffolding/test_polish/isoseq_polished/merge_anno_bed_to_fa

```

# Getting open read frames (ORF) from the transcript sequences

```{bash}

python /home/jc748673/tama/tama_go/orf_nmd_predictions/tama_orf_seeker.py -f /home/jc748673/Scaffolding/test_polish/isoseq_polished/merge_anno_bed_to_fa -o /home/jc748673/Scaffolding/test_polish/isoseq_polished/orf_aa_merged_polished

```

# Split orf_aa_merged_polished into files for blastp (https://github.com/GenomeRIK/tama/wiki/TAMA-GO:-Split-Files). The number of outputs is automatically chosen by tama_fasta_splitter or can be selected by the user. In this case, split into 9 output files. 

```{bash}

python /home/jc748673/tama/tama_go/split_files/tama_fasta_splitter.py /home/jc748673/Scaffolding/test_polish/isoseq_polished/orf_aa_merged_polished /home/jc748673/Scaffolding/test_polish/isoseq_polished/split_orf_aa_merged_polished num_splits

```
# Make the Uniprot/Uniref/Whatever protein database file to a blastp database file

```{bash}

makeblastdb -in /home/jc748673/Scaffolding/test_polish/isoseq_polished/uniprot_sprot.fasta -dbtype prot

```


# Blasting the amino acid sequences against the Uniprot/Uniref (in this case it was done with UniprotKB_SwissProt) protein database (run for each of the nine  spliced files)

```{bash}

blastp -evalue 1e-10 -num_threads 16 -db /home/jc748673/Scaffolding/test_polish/isoseq_polished/uniprot_sprot.fasta -query /home/jc748673/Scaffolding/test_polish/isoseq_polished/split_orf_aa_merged_polished_1.fa

```

# Concatanate all of the nine blastp outputs


# Parsing the Blastp output file for top hits

```{bash}

python tama_orf_blastp_parser.py -b ${blastp} -o ${outfile}
python tama_orf_blastp_parser.py -b ${blastp} -o ${outfile} -f ensembl

```

# Create new bed file with CDS regions

python tama_cds_regions_bed_add.py -p ${orf} -a ${bed} -f ${fasta} -o ${outfile}

