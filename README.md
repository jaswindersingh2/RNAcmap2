# RNAcmap2
An improved fully automatic pipeline for predicting contact maps of RNAs by evolutionary coupling analysis

## System Requirments

**Hardware Requirments:**
It is recommended that your system should have 32 GB RAM, 500 GB disk space to support the in-memory operations for RNA sequence length less than 500. Multiple CPU threads are also recommended as the MSA generating process is computationally expensive.

**Software Requirments:**
* [Python3.6](https://docs.python-guide.org/starting/install3/linux/)
* [Perl-5.4 or later](https://www.perl.org/get.html)
* [Anaconda or Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
* [gawk 4.1 or later](https://www.gnu.org/software/gawk)

RNAcmap2 has been tested on Ubuntu 14.04, 16.04, and 18.04 operating systems.

## Installation of RNAcmap and its dependencies

Clone RNAcmap2 github repo:

1. `git clone https://github.com/jaswindersingh2/RNAcmap2.git && cd RNAcmap2`

To create **Conda** virtual environment:

2. `conda create -n venv_rnacmap2 python=3.6 && conda activate venv_rnacmap2`

To install dependencies:

3. `while read p; do conda install --yes $p; done < requirements.txt`

To install **RNAfold** predictor for base pair probability features:

4. `conda install -c bioconda viennarna`

To install **BLAST-N** and **INFERNAL** tools for mulitple-sequence-alignment search:

5. `conda install -c bioconda blast`
6. `conda install -c bioconda infernal`
7. `conda install -c bioconda easel`
8. `conda install -c bioconda seqkit`

To setup `utils` scripts

9. link seqkit executable (`$conda_dir/envs/venv_rnacmap2/bin/seqkit`) to `RNAcmap2/utils/seqkit`
10. install [HHsuite-github-repo](https://github.com/soedinglab/hh-suite/tree/master/scripts) and link `$hh-suite_dir/scripts/reformat.pl` to `RNAcmap2/utils/reformat.pl`
11. install [RNAsol standalone program](https://yanglab.nankai.edu.cn/RNAsol/) and link `$RNAsol_dir/package/bin/getpssm.pl` to `RNAcmap2/utils/getpssm.pl`, link `parse_blastn_local.pl` in the same way.
12. install PLMC ([plmc-github-repo](https://github.com/debbiemarkslab/plmc)) and make sure you can find plmc executable  at `RNAcmap2/plmc/bin/plmc`


## Usage


### To run RNAcmap2:

`./run_rnacmap2.sh sample_input/1eiy_C`

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

Licence
-----
Mozilla Public License 2.0


Contact
-----
jaswinder.singh3@griffithuni.edu.au, yaoqi.zhou@griffith.edu.au
