# Iso-Seq analysis (based on https://ucdavis-bioinformatics-training.github.io/2020-september-isoseq/liz/bioconda/2-bioconda)

# Load PacBio (bioconda) tools

module load  pbbioconda/20231004 #for HPC

# Merge iso-seq flnc from all samples/tissues

ls H_liver-flnc.bam H_gonad-flnc.bam H_eyes-flnc.bam H_gills-flnc.bam H_tentacle-flnc.bam > flnc_all.fofn

isoseq3 cluster flnc_all.fofn clustered.bam --use-qvs

# Align Full Length Non Chemiric (FLNC) reads to ref genome
pbmm2 align /home/jc748673/Scaffolding/yahs/yahs.out_scaffolds_final.fa /home/jc748673/Genome/IsoSeq/IsoSeq_Analysis/clustered.hq.bam aln.all_hq.bam --sort  --preset ISOSEQ  --log-level INFO

# Collapse sorted aligned isoseq (bam) 
isoseq3 collapse aln.all.bam alz.collapsed.all.gff
