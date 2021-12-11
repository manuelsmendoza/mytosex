# Description
`MyToSex` is  a novel tool for *in silico* sex determination based on the mitochondrial genomes content. 
Some mussels and clams species have an unusual system of mitochondrial inheritance termed double uniparental 
inheritance, which involves the transmission of two different sex-associated mitogenomes haplotypes to the offspring.
Females contain only F-type mitogenomes (mtF) whereas males carry both haplotypes mtF, and also M-type (mtM). This
tools works in three acts:

1. **Mitogenomes detection and quantification**: To detect the mitogenomes presence, we mapped all the reads to both 
   mitotypes of reference. Afterward, we calculate the mitogenomes coverage, sequencing depth and depth uniformity.

2. **Sex prediction**: We implemented an artificial neural network that takes those three parameters of each mitogenome 
   as input to predict the sex of the sample.

3. **Supporting analysis**: We also implemented two additional analyses to complement and bring more support the 
   results.
   1. We clustered the samples applying dimensional reduction techniques to verify if the resultant clustering agree 
      with the prediction.
   2. We use the reads mapped to the mitogenomes to assemble *de novo* the mitogenes and after annotated them, we 
      performed a phylogenetic analysis building the different gene trees.

# Installation and use
`MyToSex` is an open-source tools written in [Python3](https://www.python.org) that requires multiple modules.
This tool can be installed as follows.

```shell
# Create the conda environment
conda env create --file environment.yaml

# Clone the repository
git clone https://github.com/manuelsmendoza/mytosex.git

# Add it to the PATH
export PATH=$PATH:$(readlink -f mytosex)
```

To use this tools, create a file containing the information required for the analysis and pass it as positional 
argument as we have specified below. To see more information about the settings file, check it in the 
[wiki](https://github.com/manuelsmendoza/mytosex/wiki).
```shell
python3 mytosex.py settings.yaml
```

# Resources
This tool was tested with the following resources:
- High-performance computer:
  - OS: Gentoo Base System release 2.7 x86_64 
  - CPU: Intel Xeon E5-2680 v3 (24) @ 3.300GHz
  - GPU: ASPEED Technology, Inc. ASPEED Graphics Family
  - Memory: 10237MiB / 128666MiB
- Personal computer:
  - OS: Ubuntu 18.04.6 LTS x86_64
  - CPU: AMD Ryzen 9 3950X (32) @ 3.500GHz
  - GPU: AMD ATI 10:00.0 Device 731f 
  - Memory: 697MiB / 32041MiB

# Citation
If you only use `MyToSex` cite us as follows:

Mendoza M. and Canchaya A., MyToSex: Sexual inference based on mitochondrial genome content [...]

Please, also include to:

- Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie2. *Nat Methods*. 2012;9(4):357-359.
    doi:[10.1038/nmeth.1923](https://www.nature.com/articles/nmeth.1923).
    ```text
    @article{langmead2012bowtie2, 
      title={Fast gapped-read alignment with Bowtie2},
      author={Langmead, Ben and Salzberg, Steven L},
      journal={Nature methods},
      volume={9},
      number={4},
      pages={357--359},
      year={2012},
      publisher={Nature Publishing Group}
    }
    ```

- Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. *Bioinformatics*. 
    2009;25(16):2078-2079. 
    doi:[10.1093/bioinformatics/btp352](https://academic.oup.com/bioinformatics/article/25/16/2078/204688).
    ```text
    @article{li2009samtools,
      title={The sequence alignment/map format and SAMtools},
      author={Li, Heng and Handsaker, Bob and Wysoker, Alec and Fennell, Tim and Ruan, Jue and Homer, Nils and Marth, Gabor and Abecasis, Goncalo and Durbin, Richard},
      journal={Bioinformatics},
      volume={25},
      number={16},
      pages={2078--2079},
      year={2009},
      publisher={Oxford University Press}
    }
    ```

If you also perform the supporting analysis, please cite them too.
- Samples clustering:
  - McInnes, L., Healy, J. and Melville, J. UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction.
    *arXiv preprint*, arXiv:1802.03426.
    ```text
    @misc{mcinnes2020umap,
      title={UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction}, 
      author={Leland McInnes and John Healy and James Melville},
      year={2020},
      eprint={1802.03426},
      archivePrefix={arXiv}
    } 
    ```

- Phylogenetic analysis:
  - Haas BJ, Papanicolaou A, Yassour M, *et. al*. *De novo* transcript sequence reconstruction from RNA-seq using the 
    Trinity platform for reference generation and analysis. *Nat Protoc*. 2013;8(8):1494-1512. 
    doi:[10.1038/nprot.2013.084](https://www.nature.com/articles/nprot.2013.084) 
    ```text
    @article{haas2013trinity,
      title={De novo transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis},
      author={Haas, Brian J and Papanicolaou, Alexie and Yassour, Moran and Grabherr, Manfred and Blood, Philip D and Bowden, Joshua and Couger, Matthew Brian and Eccles, David and Li, Bo and Lieber, Matthias and others},
      journal={Nature protocols},
      volume={8},
      number={8},
      pages={1494--1512},
      year={2013},
      publisher={Nature Publishing Group}
    }
    ```
  - Camacho C, Coulouris G, Avagyan V, *et al*. BLAST+: architecture and applications. *BMC Bioinformatics*. 
    2009;10:421. 
    doi:[10.1186/1471-2105-10-421](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-421).
    ```text
    @article{camacho2009blast+,
      title={BLAST+: architecture and applications},
      author={Camacho, Christiam and Coulouris, George and Avagyan, Vahram and Ma, Ning and Papadopoulos, Jason and Bealer, Kevin and Madden, Thomas L},
      journal={BMC bioinformatics},
      volume={10},
      number={1},
      pages={421--429},
      year={2009},
      publisher={Springer}
    }
    ```
  - Katoh K, Kuma K, Toh H, Miyata T. MAFFT version 5: improvement in accuracy of multiple sequence alignment. 
    *Nucleic Acids Res*. 2005;33(2):511-518. 
    doi:[10.1093/nar/gki198](https://academic.oup.com/nar/article/33/2/511/2549118).
    ```text
    @article{katoh2005mafft,
      title={MAFFT version 5: improvement in accuracy of multiple sequence alignment},
      author={Katoh, Kazutaka and Kuma, Kei-ichi and Toh, Hiroyuki and Miyata, Takashi},
      journal={Nucleic acids research},
      volume={33},
      number={2},
      pages={511--518},
      year={2005},
      publisher={Oxford University Press}
    }
    ```
  - Darriba D, Posada D, Kozlov AM, Stamatakis A, Morel B, Flouri T. ModelTest-NG: A New and Scalable Tool for the 
    Selection of DNA and Protein Evolutionary Models. *Mol Biol Evol*. 2020;37(1):291-294. 
    doi:[10.1093/molbev/msz189](https://academic.oup.com/mbe/article/37/1/291/5552155).
    ```text
    @article{darriba2020modeltest,
      title={ModelTest-NG: a new and scalable tool for the selection of DNA and protein evolutionary models},
      author={Darriba, Diego and Posada, David and Kozlov, Alexey M and Stamatakis, Alexandros and Morel, Benoit and Flouri, Tomas},
      journal={Molecular Biology and Evolution},
      volume={37},
      number={1},
      pages={291--294},
      year={2020},
      publisher={Oxford University Press}
    }
    ```
  - Stamatakis A. RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. 
    *Bioinformatics*. 2014;30(9):1312-1313. 
    doi:[10.1093/bioinformatics/btu033](https://academic.oup.com/bioinformatics/article/30/9/1312/238053).
    ```text
    @article{stamatakis2014raxml,
      title={RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies},
      author={Stamatakis, Alexandros},
      journal={Bioinformatics},
      volume={30},
      number={9},
      pages={1312--1313},
      year={2014},
      publisher={Oxford University Press}
    }
    ```

# Acknowledgement
This work was supported by the European Social Fund and the Government of Xunta de Galicia (Scholarship reference 
ED481A-2018/305 awarded by Manuel Mendoza).

We developed this tools using the computational resources of the 
[Supercomputing Center of Galicia (CESGA)](https://www.cesga.es) using [Pycharm](https://www.jetbrains.com/pycharm/) 
with an Academic License freely provided for JetBrain.
