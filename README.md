# Description
`MyToSex` is  a novel tool for *in silico* sex determination based on the mitochondrial genomes content. 
Some mussels and clams species have an unusual system of mitochondrial inheritance termed double uniparental 
inheritance, which involves the transmission of two different sex-associated mitogenomes haplotypes to the offspring.
Females contain only F-type mitogenomes (mtF) whereas males carry both haplotypes mtF, and also M-type (mtM). This
tools works in two different acts:

1. **Mitogenomes detection and quantification**: To detect the mitogenomes presence, we mapped all the reads to both
    mitotypes. From this alignment we extracted some metrics that are used to determine the sex.

2. **Additional analyses**: We also implemented two additional analyses to complement and bring more support the 
   results.
   1. Samples clustering: We extracted multiples metrics from the reads alignments to the mitogenomes and applying 
      a dimensional reduction ([UMAP](https://arxiv.org/abs/1802.03426)) to verify if the resultant clustering agree 
      with the sex-determination results obtained previously.
   2. Phylogenetic analysis of protein-coding mitogenes: We use the reads that mapped to the mitogenomes to assemble
      *de novo* the protein-coding genes which are used to perform a phylogenetic analysis incorporating the mitogenes
      of reference and also adding information from other species.

# Installation
`MyToSex` is an open-source tools written in Python3 that requires the following modules: ... Furthermore, it also calls
third-party software.

# Citation

# Acknowledgement
This work was supported by the European Social Fund and the Government of Xunta de Galicia (Scholarship reference 
ED481A-2018/305 awarded by Manuel Mendoza).

We developed this tools using the computational resources of the 
[Supercomputing Center of Galicia (CESGA)](https://www.cesga.es) using [Pycharm](https://www.jetbrains.com/pycharm/) 
with an Academic License freely provided for JetBrain.

