# Assembly of mitochondrial genome using Pacbio HiFi reads and/or Hifiasm assembly with [MitoHiFi](https://github.com/marcelauliano/MitoHiFi)

**MitoHiFi** allows the assembly of mitochondrial genomes with raw Pacbio HiFi reads (flag -r) or from the assembled contigs/scaffolds (flag -c).
In addition, you should download the mitochondrial genome assembly of the most closely related species (fasta file and genebank file are required). This can be done by using findMitoReference.py script (provided with the tool). 
There is a great YouTube video (**BGA 2023**) by the author of the tool (**Marcela Uliano-Silva**) which is recommended to watch: 
[MitoHiFi_BGA2023](https://www.youtube.com/watch?v=1NWHC2zkRmg&t=1036s)

## This is the line of code used for the abalone mitochondrial genome assembly

You should change/add parameters according to your available data and type of organism. For example, for abalone, it was required to use the genetic code (-o) for invertebrates (5). 
For other specifications look at https://github.com/marcelauliano/MitoHiFi

```
mitohifi.py -c scaffolds_final.fa -f ref_mito_genome.fasta -g ref_mito_genome.gb -t 16 -a animal -o 5
```
The output includes multiple files including stats, a list of genes, a mito genome plot and more. 

## Genetic codes list (use the most appropriate code according to your organism) : 

1. The Standard Code **-o 1** 
2.  The Vertebrate Mitochondrial Code **-o 2**
3.  The Yeast Mitochondrial Code **-o 3**
4.  The Mold,Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code **-o 4**
5.  The Invertebrate Mitochondrial Code **-o 5**
6.  The Ciliate, Dasycladacean and Hexamita Nuclear Code **-o 6**
7.  The Echinoderm and Flatworm Mitochondrial Code **-o 9**
8.  The Euplotid Nuclear Code **-o 10**
9.  The Bacterial, Archaeal and Plant Plastid Code **-o 11**
10. The Alternative Yeast Nuclear Code **-o 12**
11. The Ascidian Mitochondrial Code **-o 13**
12. The Alternative Flatworm Mitochondrial Code **-o 14**
13. The Chlorophycean Mitochondrial Code **-o 16**
14. The Trematode Mitochondrial Code **-o 21** 
15. The Scenedesmus obliquus Mitochondrial Code **-o 22**
16. The Thraustochytrium Mitochondrial Code **-o 23**
17. The Pterobranchia Mitochondrial Code **-o 24**
18. The Candidate Division SR1 and Gracilibacteria Code **-o 25**
