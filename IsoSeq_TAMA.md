# map align isoseq high-quality isoforms to genome

```{bash}

minimap2 -ax splice -uf -C5  /home/jc748673/Scaffolding/test_polish/yahs_polished.out_scaffolds_final.fa /home/jc748673/Genome/IsoSeq/IsoSeq_Analysis/High_Quality_Isoforms.fasta > aln_iso_hq_polished.sam

```

# sort sam

samtools sort aln_iso_hq_polished.sam -O sam -o aln_iso_hq_polished_sort.sam

# https://github.com/GenomeRIK/tama/wiki/TAMA-GO:-ORF-and-NMD-predictions (link to github manual TAMA)

# tama collapse capped

singularity run /fast/tmp/containers/python-2.7.17.sif python /home/jc748673/tama/tama_collapse.py -s /home/jc748673/Scaffolding/test_polish/isoseq_polished/aln_iso_hq_polished_sort.sam -f /home/jc748673/Scaffolding/test_polish/yahs_polished.out_scaffolds_final.fa -p /home/jc748673/Scaffolding/test_polish/isoseq_polished/iso_hq_tama_polished -x no_cap

# tama collapse capped

singularity run /fast/tmp/containers/python-2.7.17.sif python /home/jc748673/tama/tama_collapse.py -s /home/jc748673/Scaffolding/test_polish/isoseq_polished/aln_iso_hq_polished_sort.sam -f /home/jc748673/Scaffolding/test_polish/yahs_polished.out_scaffolds_final.fa -p /home/jc748673/Scaffolding/test_polish/isoseq_polished/iso_hq_tama_polished_capped -x capped

# tama_degradation_signature

singularity run /fast/tmp/containers/python-2.7.17.sif python /home/jc748673/tama/tama_go/file_stats/tama_degradation_signature.py -c /home/jc748673/Scaffolding/test_polish/isoseq_polished/iso_hq_tama_polished_capped_trans_read.bed -nc /home/jc748673/Scaffolding/test_polish/isoseq_polished/iso_hq_tama_polished_trans_read.bed -o /home/jc748673/Scaffolding/test_polish/isoseq_polished/tama_degradation_signature

# tama merge bed files annotations

python /home/jc748673/tama/tama_merge.py -f /home/jc748673/Scaffolding/test_polish/isoseq_polished/filelist.txt -p /home/jc748673/Scaffolding/test_polish/isoseq_polished/merged_annos

# filelist.txt for merging beds:

/home/jc748673/Scaffolding/test_polish/isoseq_polished/iso_hq_tama_polished_capped.bed	capped	1,1,1	caplib
/home/jc748673/Scaffolding/test_polish/isoseq_polished/iso_hq_tama_polished.bed	capped	2,1,1	nocaplib

# Converting the bed files into fasta files

bedtools getfasta -name -split -s -fi /home/jc748673/Scaffolding/test_polish/yahs_polished.out_scaffolds_final.fa -bed /home/jc748673/Scaffolding/test_polish/isoseq_polished/merged_annos.bed -fo /home/jc748673/Scaffolding/test_polish/isoseq_polished/merge_anno_bed_to_fa

# Getting open read frames (ORF) from the transcript sequences

python /home/jc748673/tama/tama_go/orf_nmd_predictions/tama_orf_seeker.py -f /home/jc748673/Scaffolding/test_polish/isoseq_polished/merge_anno_bed_to_fa -o /home/jc748673/Scaffolding/test_polish/isoseq_polished/orf_aa_merged_polished

# Blasting the amino acid sequences against the Uniprot/Uniref protein database

blastp -evalue 1e-10 -num_threads 16 -db https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100 -query /home/jc748673/Scaffolding/test_polish/isoseq_polished/orf_aa_merged_polished

# Parsing the Blastp output file for top hits

python tama_orf_blastp_parser.py -b ${blastp} -o ${outfile}
python tama_orf_blastp_parser.py -b ${blastp} -o ${outfile} -f ensembl

# Create new bed file with CDS regions

python tama_cds_regions_bed_add.py -p ${orf} -a ${bed} -f ${fasta} -o ${outfile}

