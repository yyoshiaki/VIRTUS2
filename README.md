# VIRTUS2 : VIRal Transcript Usage Sensor v2.0 <img src="https://github.com/yyoshiaki/VIRTUS/raw/master/img/VIRTUS.jpg" width="20%" align="right" />

**Note : We updated VIRTUS to version2. In this version, we removed the gene quantification step by Salmon, and single virus mode to focus on virus-wide exploration.**

Virus transcript detection and quantification using normal human RNAseq. VIRTUS is the first tool to detect viral transcripts considering their splicing event rather than the viral genome copy number. VIRTUS can be applied to both bulk RNAseq and single-cell RNAseq. The virus reference covers 762 viruses including SARS-CoV-2 (cause of COVID-19). The workflow is implemented by [Common Workflow Language](https://www.commonwl.org/) and [Rabix](https://rabix.io/). You can specify each parameter individually or give `yaml` or `json` file which describes all the parameter information. In detail, check [the CWL User Guide](http://www.commonwl.org/user_guide/) out. 

![img](https://github.com/yyoshiaki/VIRTUS/raw/master/img/webimage.jpg)

## Contact

Yoshiaki Yasumizu ([yyasumizu@ifrec.osaka-u.ac.jp](yyasumizu@ifrec.osaka-u.ac.jp))

## Citation

VIRTUS: a pipeline for comprehensive virus analysis from conventional RNA-seq data
Yasumizu, Yoshiaki, Atsushi Hara, Shimon Sakaguchi, and Naganari Ohkura. 2020. “OUP Accepted Manuscript.” Edited by Jan Gorodkin. *Bioinformatics*, October. https://doi.org/10.1093/bioinformatics/btaa859.

## Acknowledgement

- [Kozo Nishida](https://github.com/kozo2) : CI config and many other supports.
- Ayano Onishi : Logo design.

## Licence

This software is freely available for academic users. Usage for commercial purposes is not allowed. Please refer to the LICENCE page.

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a>

## Tutorial


### Installation

1. Install cwltool and python packages.

If you are not using python3, install python3 via anaconda or another way. Then, install these;

```
pip install cwlref-runner numpy pandas scipy statsmodels seaborn
```

2. Install Docker

Install docker (alternatively, you can use udocker when you have no root privileges. See [Tips section](https://github.com/yyoshiaki/VIRTUS2#tips).)

3. clone the repo and add the path.

```
git clone https://github.com/yyoshiaki/VIRTUS2
```

You can add `./VIRTUS/bin` to `PATH` in `.zshrc` or `.bashrc` etc. 


### Create indices

just run the code below. A set of indices will be created in your current directory. You can specify a yaml file, options or just a url as the input.

```
createindex.cwl https://raw.githubusercontent.com/yyoshiaki/VIRTUS2/master/workflow/createindex.job.yaml
```

### First run

VIRTUS has a great wrapper for multiple samples. The input is a comma-separated text file or CSV file. The first column is arbitral sample names, the second is SRR id, or fastq files (when you specify `--fastq` option). Note that `--fastq` requires the suffix removed file names. Refer to the documentation in more detail.  The third column is the sequence layout (SE or PE), and the Fourth is groups. Let's create an example csv file (or [download it](https://raw.githubusercontent.com/yyoshiaki/VIRTUS2/master/wrapper/input.csv)).

input.csv
```
Name,SRR,Layout,Group
Flu_1,SRR9856912,PE,H3N2
Flu_2,SRR9856913,PE,H3N2
Ctrl_1,SRR9856914,PE,Mock
Ctrl_2,SRR9856915,PE,Mock
```

Then, run this (edit DIR_INDEX_ROOT).

```
DIR_INDEX_ROOT=/dir/to/indeices/created/
VIRTUS_wrapper.py input.csv \
    --genomeDir_human $DIR_INDEX_ROOT/STAR_index_human \
    --genomeDir_virus $DIR_INDEX_ROOT/STAR_index_virus \
    --nthreads=4
```

You can get this heatmap and `summary.csv` which contains the ratio of viral reads (hit viral reads/read mapped on the human genome), the stats by Mann-Whitney's U-test, and its false discovery rate.

![img](https://github.com/yyoshiaki/VIRTUS2/blob/master/img/clustermap.png)

## Tips