# RNAcmap2
An improved fully automatic pipeline for predicting contact maps of RNAs by evolutionary coupling analysis


## System Requirments

**Hardware Requirments:**
It is recommended that your system should have 32 GB RAM, 500 GB disk space to support the in-memory operations for RNA sequence length less than 500. Multiple CPU threads are also recommended as the MSA generating process is computationally expensive.

**Software Requirments:**
* [Python3.6](https://docs.python-guide.org/starting/install3/linux/)
* [Perl-5.4 or later](https://www.perl.org/get.html)
* [Anaconda or Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

RNAcmap2 has been tested on Ubuntu 14.04, 16.04, and 18.04 operating systems.


## Installation of RNAcmap2 and its dependencies

Clone RNAcmap2 github repo:

1. `git clone https://github.com/jaswindersingh2/RNAcmap2.git && cd RNAcmap2`

Just run the following command to create **Conda** virtual environment and install **Conda** dependencies:

2. `conda env create --file environment.yaml`

3. `conda activate venv_rnacmap2`


## Install DCA predictor using following commands:

For mfDCA and plmDCA:

4. `pip install pydca`

For PLMC:

5. `git clone https://github.com/debbiemarkslab/plmc && cd plmc && make all-openmp && cd -`

For GREMLIN:

6. `git clone "https://github.com/sokrypton/GREMLIN_CPP" && cd GREMLIN_CPP && g++ -O3 -std=c++0x -o gremlin_cpp gremlin_cpp.cpp -fopenmp && cd ../`


## Download the reference database used by RNAcmap2 using following command:

7. `./db_download.sh`


## To format the database to use with **blastn**, the following command can be used.

8. `makeblastdb -in ./database/nt_metagenomics_database/nt_metagenomics2 -dbtype nucl`


## Usage


### To run RNAcmap2:

9. `./run_rnacmap2.sh 6p2h_A mfdca`

## Reproduce results of RNAcmap pipeline:

Refer to [benchmarking](https://github.com/jaswindersingh2/RNAcmap2/tree/main/benchmarking) folder of this repo.

## Third party programs

* cmbuild, cmcalibrate, and cmsearch from [INFERNAL tool](http://eddylab.org/infernal) version 1.1.4
* esl-reformat from [easel tool](https://anaconda.org/bioconda/easel) version 0.48
* blastn and makeblastdb from [BLAST tool](https://anaconda.org/bioconda/blast) version 2.11.0
* RNAfold from [ViennaRNA](https://anaconda.org/bioconda/viennarna) version 2.4.18
* utils/reformat.pl from [HHsuite-github-repo](https://github.com/soedinglab/hh-suite/tree/master/scripts)
* utils/getpssm.pl and utils/parse\_blastn\_local.pl from [RNAsol standalone program](https://yanglab.nankai.edu.cn/RNAsol/)
* utils/seqkit from [seqkit toolkit](https://bioinf.shenwei.me/seqkit/)
* PLMC from [plmc-github-repo](https://github.com/debbiemarkslab/plmc)
* GREMLIN from [gremlin-github-repo](https://github.com/sokrypton/GREMLIN_CPP)
* mfDCA and plmDCA from [pydca-github-repo](https://github.com/KIT-MBS/pydca)


## Citation guide

**If use RNAcmap2 for your research, please cite the following papers:**

Jaswinder Singh, Thomas Litfin, Kuldip Paliwal, Jaspreet Singh, and Yaoqi Zhou. "An improved fully automatic pipeline for predicting contact maps of RNAs by evolutionary coupling analysis."

**If use RNAcmap2 pipeline, please consider citing the following papers:**

BLAST-N:

[1] Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J., Zhang, Z., Miller, W. and Lipman, D.J., 1997. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucleic acids research, 25(17), pp.3389-3402.


INFERNAL:

[2] Nawrocki, E.P. and Eddy, S.R., 2013. Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics, 29(22), pp.2933-2935.


RNAfold:

[3] Lorenz, R., Bernhart, S.H., Zu Siederdissen, C.H., Tafer, H., Flamm, C., Stadler, P.F. and Hofacker, I.L., 2011. ViennaRNA Package 2.0. Algorithms for molecular biology, 6(1), pp.1-14.


RNAcmap Pipeline:

[4] Zhang, T., Singh, J., Litfin, T., Zhan, J., Paliwal, K. and Zhou, Y., 2021. RNAcmap: a fully automatic pipeline for predicting contact maps of RNAs by evolutionary coupling analysis. Bioinformatics.


PLMC:

[5] Hopf, T.A., Ingraham, J.B., Poelwijk, F.J., Schärfe, C.P., Springer, M., Sander, C. and Marks, D.S., 2017. Mutation effects predicted from sequence co-variation. Nature biotechnology, 35(2), pp.128-135.


GREMLIN:

[6] Kamisetty, H., Ovchinnikov, S. and Baker, D., 2013. Assessing the utility of coevolution-based residue–residue contact predictions in a sequence-and structure-rich era. Proceedings of the National Academy of Sciences, 110(39), pp.15674-15679.


mfDCA and plmDCA:

[7] Zerihun, MB., Pucci, F, Peter, EK, and Schug, A. pydca: v1.0: a comprehensive software for direct coupling analysis of RNA and protein sequences. Bioinformatics, btz892, doi.org/10.1093/bioinformatics/btz892

[8] Morcos, F., Pagnani, A., Lunt, B., Bertolino, A., Marks, DS., Sander, C., Zecchina, R., Onuchic, JN., Hwa, T., and Weigt, M. Direct-coupling analysis of residue coevolution captures native contacts across many protein families PNAS December 6, 2011 108 (49) E1293-E1301, doi:10.1073/pnas.1111471108

[9] Ekeberg, M., Lövkvist, C., Lan, Y., Weigt, M., & Aurell, E. (2013). Improved contact prediction in proteins: Using pseudolikelihoods to infer Potts models. Physical Review E, 87(1), 012707, doi:10.1103/PhysRevE.87.012707


SeqKit:

[10] Shen, W., Le, S., Li, Y. and Hu, F., 2016. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PloS one, 11(10), p.e0163962.


**If use RNAcmap2 datasets, please consider citing the following papers:**

Protein Data Bank (PDB):

[11] Berman, H.M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T.N., Weissig, H., Shindyalov, I.N. and Bourne, P.E., 2000. The protein data bank. Nucleic acids research, 28(1), pp.235-242.

CD-HIT-EST:

[12] Fu, L., Niu, B., Zhu, Z., Wu, S. and Li, W., 2012. CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics, 28(23), pp.3150-3152.


Licence
-----
Mozilla Public License 2.0


Contact
-----
jaswinder.singh3@griffithuni.edu.au, yaoqi.zhou@griffith.edu.au
