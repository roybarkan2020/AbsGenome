# Single-sample telomere-to-telomere assembly with HiFi, ultralong and Hi-C reads

singularity run /sw/containers/hifiasm-0.18.2r467beta.sif hifiasm -o Hal_Asi.asm --h1 Omni1.fq.gz --h2 Omni2.fq.gz --ul ul.fq.gz Hal_Asi_HiFi.fq.gz

# gfa to fasta (hifiasm primery assembly) 
awk '/^S/{print ">"$2;print $3}' Hal_Asi.asm.bp.p_ctg.gfa   > Hal_Asi.asm.bp.p_ctg.fa
