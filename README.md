# VIRTUS2 : VIRal Transcript Usage Sensor v2.0.2 <img src="https://github.com/yyoshiaki/VIRTUS/raw/master/img/VIRTUS.jpg" width="20%" align="right" />

**!!Note : We updated VIRTUS to version2. In this version, we removed the gene quantification step by Salmon and single virus mode, and added coverage on viral genomes to the result to focus on virus-wide exploration. If you want to use the single virus mode, visit [https://github.com/yyoshiaki/VIRTUS](https://github.com/yyoshiaki/VIRTUS)**

Virus transcript detection and quantification using normal human RNAseq. VIRTUS is the first tool to detect viral transcripts considering their splicing event rather than the viral genome copy number. VIRTUS can be applied to both bulk RNAseq and single-cell RNAseq. The virus reference covers 762 viruses including SARS-CoV-2 (cause of COVID-19). The workflow is implemented by [Common Workflow Language](https://www.commonwl.org/) and [Rabix](https://rabix.io/). You can specify each parameter individually or give `yaml` or `json` file which describes all the parameter information. In detail, check [the CWL User Guide](http://www.commonwl.org/user_guide/) out. 

![img](https://github.com/yyoshiaki/VIRTUS2/raw/master/img/webimage.jpg)

## Update from [VIRTUS1](https://github.com/yyoshiaki/VIRTUS)

To focus on virus-wide exploration, we changed several functions from VIRTUS1.

- removed the single virus mode
- added coverage on the viral genome and quality information on the results.
- removed the default filter for low hit viruses.
- new visualization by VIRTUS_wrapper.py

## Tutorials

- [Install](#install)
- [Tools](#tools)
- [Tips](#tips)
- [VIRTUS2 for single cell RNAseq](#virtus2-for-single-cell-rnaseq)

## Install

### dependencies

- python3
- cwltool `pip install cwltool` 
- docker (alternatively, you can use udocker, or singularity mode when you have no root privileges. See [Tips section](https://github.com/yyoshiaki/VIRTUS#tips).)
- STAR requires about 37Gb RAM.

### VIRTUS2

1. install cwltool

```
pip install cwltool
```

In detail, check out [the cwltool website](https://github.com/common-workflow-language/cwltool).
(Tested version is cwltool==3.1.20210922203925 and 3.1.20221008225030)

2. Setup docker

Setup [docker](https://www.docker.com/). If you are using Mac, increase the memory limit above 40Gb([Documentation](https://docs.docker.com/docker-for-mac/#resources)).

3. Clone VIRTUS2.
```
git clone https://github.com/yyoshiaki/VIRTUS2
```

You can add `./VIRTUS/bin` to `PATH` in `.zshrc` or `.bashrc` etc. If you have used VIRTUS1, you need to remove the path to VIRTUS1 to avoid conflicts.

4. Version confirmation

`Tool --help` will tell you the version. ex. `VIRTUS.PE.cwl --help`

### VIRTUS wrapper

VIRTUS wrapper uses Python3 and several packages. To install python packages, you may need to execute the command below after installing python3.

```
pip install numpy pandas scipy statsmodels seaborn
```


## Usage

1. [create indices](#create-indices)
2. run [VIRTUS2](#virtuspecwl)

## Tools

### Create indices

just run the code below. A set of indices will be created in your current directory. You can specify a yaml file, options or just a url as the input.

```
createindex.cwl https://raw.githubusercontent.com/yyoshiaki/VIRTUS2/master/workflow/createindex.job.yaml
```

#### createindex.cwl (execute only once)

Fetch reference data and create indexes for VIRTUS2. Located in `VIRTUS2/workflow`.

```
usage: ./createindex.cwl [-h] --url_virus URL_VIRUS \ 
                              --output_name_virus OUTPUT_NAME_VIRUS \
                              [--runThreadN RUNTHREADN] \
                              --dir_name_STAR_virus DIR_NAME_STAR_VIRUS \
                              --url_genomefasta_human URL_GENOMEFASTA_HUMAN \
                              --output_name_genomefasta_human OUTPUT_NAME_GENOMEFASTA_HUMAN \ 
                              --dir_name_STAR_human DIR_NAME_STAR_HUMAN \
[job_order]

VIRTUS v2.0

positional arguments:
  job_order             Job input json file

optional arguments:
  -h, --help            show this help message and exit
  --url_virus URL_VIRUS
  --output_name_virus OUTPUT_NAME_VIRUS
  --runThreadN RUNTHREADN
  --dir_name_STAR_virus DIR_NAME_STAR_VIRUS
  --url_genomefasta_human URL_GENOMEFASTA_HUMAN
  --output_name_genomefasta_human OUTPUT_NAME_GENOMEFASTA_HUMAN
  --dir_name_STAR_human DIR_NAME_STAR_HUMAN
```

```
./createindex.cwl createindex.job.yaml
```

virus fasta is from [VirTect](https://github.com/WGLab/VirTect).

![img/createindex.png](https://github.com/yyoshiaki/VIRTUS2/raw/master/img/createindex.png)

### VIRTUS.PE.cwl

The main VIRTUS pipeline for paired-end RNA-seq. Located in `VIRTUS2/workflow`.

```
usage: ./VIRTUS.PE.cwl [-h] \
                        --fastq2 FASTQ2 \
                        --fastq1 FASTQ1 \
                        --genomeDir_human \
                        GENOMEDIR_HUMAN \
                        [--outFileNamePrefix_human OUTFILENAMEPREFIX_HUMAN] \
                        [--nthreads NTHREADS] \
                        --genomeDir_virus GENOMEDIR_VIRUS \
                        [--kz_threshold KZ_THRESHOLD] \
                        [--filename_output FILENAME_OUTPUT] [job_order]

VIRTUS v2.0

positional arguments:
  job_order             Job input json file

optional arguments:
  -h, --help            show this help message and exit
  --fastq2 FASTQ2
  --fastq1 FASTQ1
  --genomeDir_human GENOMEDIR_HUMAN
  --outFileNamePrefix_human OUTFILENAMEPREFIX_HUMAN
  --nthreads NTHREADS
  --genomeDir_virus GENOMEDIR_VIRUS
  --kz_threshold KZ_THRESHOLD
                        default : 0.1
  --filename_output FILENAME_OUTPUT
                        default : VIRTUS.output.tsv
```

kz_threshold is the parameter for the filteration of low complexity sequences. refer to [https://github.com/eclarke/komplexity](https://github.com/eclarke/komplexity) for the detail.

example1

```
./VIRTUS.PE.cwl VIRTUS.PE.job.yaml
```

example2

```
./VIRTUS.PE.cwl \
--fastq1 ../test/ERR3240275/ERR3240275_1.fastq.gz \
--fastq2 ../test/ERR3240275/ERR3240275_2.fastq.gz \
--genomeDir_human ../test/STAR_index_human \
--genomeDir_virus ../test/STAR_index_virus \
--outFileNamePrefix_human human \
--nthreads 40
```

#### Output

 The example of the output is like below.

|virus|startpos|endpos|numreads|covbases|coverage|meandepth|meanbaseq|meanmapq|rate_hit|
|--|--|--|--|--|--|--|--|--|--|
|NC_007605.1_Human_herpesvirus_4_complete_wild_type_genome|1|171823|3545.5|31323|18.2298|4.0882|38.8|194.0|0.0005228124196814582|
|NC_009334.1_Human_herpesvirus_4,_complete_genome|1|172764|278.5|8872|5.135330000000001|0.319227|38.9|2.95|4.1067059337550734e-05|
|NC_003977.1_Hepatitis_B_virus,_complete_genome|1|3215|19.0|187|5.81649|0.888336|40.1|255.0|2.8017024323643226e-06|
|NC_001716.2_Human_herpesvirus_7,_complete_genome|1|153080|2.0|317|0.207081|0.00228639|32.0|1.5|2.9491604551203395e-07|
|NC_001669.1_Simian_virus_40,_complete_genome|1|5243|2.0|234|4.46309|0.07171469999999999|37.3|255.0|2.9491604551203395e-07|
|NC_001348.1_Human_herpesvirus_3,_complete_genome|1|124884|1.0|79|0.0632587|0.00126517|38.6|255.0|1.4745802275601698e-07|
|NC_001819.1_Rauscher_murine_leukemia_virus,_complete_genome|1|8282|1.0|118|1.42478|0.024028|40.3|3.0|1.4745802275601698e-07|
|NC_006998.1_Vaccinia_virus,_complete_genome|1|194711|1.0|72|0.0369779|0.0007395580000000001|34.9|255.0|1.4745802275601698e-07|
|NC_001806.1_Human_herpesvirus_1,_complete_genome|1|152261|1.0|50|0.0328383|0.0006567669999999999|39.0|255.0|1.4745802275601698e-07|

Most columns are correspinding to the output of [samtools-coverage](http://www.htslib.org/doc/samtools-coverage.html) with following modifications.

- numreads : mapped reads on each virus (in PE results, diveded by 2).
- ratio hit : reads mapped on viral genome / read mapped on human genome

Unmapped viruses are removed from the output.

#### Tool overview

![img/VIRTUS.PE.png](https://github.com/yyoshiaki/VIRTUS2/raw/master/img/VIRTUS.PE.png)

### VIRTUS.SE.cwl

The main VIRTUS pipeline for single-end RNA-seq. Located in `VIRTUS/workflow`.

```
usage: ./VIRTUS.SE.cwl [-h] \
                        --genomeDir_human GENOMEDIR_HUMAN \
                        [--outFileNamePrefix_human OUTFILENAMEPREFIX_HUMAN] \
                        [--nthreads NTHREADS] \
                        --fastq FASTQ \
                        --genomeDir_virus GENOMEDIR_VIRUS \
                        [--kz_threshold KZ_THRESHOLD] [--filename_output FILENAME_OUTPUT] [job_order]

VIRTUS v2.0

positional arguments:
  job_order             Job input json file

optional arguments:
  -h, --help            show this help message and exit
  --genomeDir_human GENOMEDIR_HUMAN
  --outFileNamePrefix_human OUTFILENAMEPREFIX_HUMAN
  --nthreads NTHREADS
  --fastq FASTQ
  --genomeDir_virus GENOMEDIR_VIRUS
  --kz_threshold KZ_THRESHOLD
                        default : 0.1
  --filename_output FILENAME_OUTPUT
                        default : VIRTUS.output.tsv
```

example1

```
./VIRTUS.SE.cwl VIRTUS.SE.job.yaml
```

example2

```
./VIRTUS.SE.cwl \
--fastq ../test/SRR8315715_2.fastq.gz \
--genomeDir_human ../test/STAR_index_human \
--genomeDir_virus ../test/STAR_index_virus \
--outFileNamePrefix_human human \
--nthreads 40
```

![img/VIRTUS.SE.png](https://github.com/yyoshiaki/VIRTUS2/raw/master/img/VIRTUS.SE.png)

### Wrapper for multiple analysis

THe wrapper script is deposited in `VIRTUS2/wrapper`.

VIRTUS2 has a wrapper for multiple samples. The input is a comma-separated text file or CSV file. The first column is arbitral sample names, the second is SRR id, or fastq files (when you specify `--fastq` option). Note that `--fastq` requires the suffix removed file names. Refer to the documentation in more detail.  The third column is the sequence layout (SE or PE), and the Fourth is groups. Let's create an example csv file (or [download it](https://raw.githubusercontent.com/yyoshiaki/VIRTUS2/master/wrapper/input.csv)).

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
DIR_INDEX_ROOT=/absolute/path/to/dir/of/indeices/created/
VIRTUS_wrapper.py input.csv \
    --genomeDir_human $DIR_INDEX_ROOT/STAR_index_human \
    --genomeDir_virus $DIR_INDEX_ROOT/STAR_index_virus \
    --nthreads=4
```

You can get this heatmap and `summary.csv` which contains the coverage on viral genomes (_cov) the ratio of viral reads (hit viral reads/read mapped on the human genome; _rate), the stats by Mann-Whitney's U-test, and its false discovery rate.

![img](https://github.com/yyoshiaki/VIRTUS2/blob/master/img/scattermap.png)
The value is the ratio of viral reads (hit viral reads/read mapped on the human genome).

#### **input**
- experiment matrix should be separated by commas (csv format).
- Only 2 groups can be tested.

**SRR mode**

|  name  |  SRR |  Layout  | Group | ... |
| ---- | ---- | - | - | - |
|  Inf_1  | SRR9856913 | PE | infected | ...|
|  Ctrl_1  |  SRR9856914  | PE  | Mock | ... |

**fastq mode**

|  name  |  fastq |  Layout  | Group | ... |
| ---- | ---- | - | - | - |
|  Inf_1  | SRR9856913 | PE | infected | ...|
|  Ctrl_1  |  SRR9856914  | PE  | Mock | ... |

- If you want to use your own fastq, add `---fastq` option. This wrapper supports only `.fastq` and `.fastq.gz`.

- fastq file specifies path excluding `.fastq.gz` or `_1.fastq.gz` and `_2.fastq.gz`. For example, `hoge/SRR1234567.fastq.gz` is described as `hoge/SRR1234567`.

- If suffix is not `.fastq.gz` or `_1.fastq.gz` and `_2.fastq.gz`, add `-s` or `-s1` and `-s2` options.

```
usage: VIRTUS_wrapper.py [-h] \
                        [--VIRTUSDir VIRTUSDIR] \
                        --genomeDir_human GENOMEDIR_HUMAN \
                        --genomeDir_virus GENOMEDIR_VIRUS \
                        [--outFileNamePrefix_human OUTFILENAMEPREFIX_HUMAN] \
                        [--nthreads NTHREADS] \
                        [-s SUFFIX_SE] \
                        [-s1 SUFFIX_PE_1] \
                        [-s2 SUFFIX_PE_2] \
                        [--fastq] \
                        [--figsize FIGSIZE] \
                        [--th_cov TH_COV] [--th_rate TH_RATE] \
                        [--singularity] \
                        input_path

positional arguments:
  input_path

optional arguments:
  -h, --help            show this help message and exit
  --VIRTUSDir VIRTUSDIR
  --genomeDir_human GENOMEDIR_HUMAN
  --genomeDir_virus GENOMEDIR_VIRUS
  --outFileNamePrefix_human OUTFILENAMEPREFIX_HUMAN
  --nthreads NTHREADS
  -s SUFFIX_SE, --Suffix_SE SUFFIX_SE
  -s1 SUFFIX_PE_1, --Suffix_PE_1 SUFFIX_PE_1
  -s2 SUFFIX_PE_2, --Suffix_PE_2 SUFFIX_PE_2
  --fastq
  --figsize FIGSIZE     (default:6,6)
  --th_cov TH_COV       threshold of max viral coverage to plot, test (default:10)
  --th_rate TH_RATE     threshold of max rate virus/human to plot, test (default:0.0001)
  --singularity         run with singularity (default:False)
```
example
```
./VIRTUS_wrapper.py input.csv \
    --VIRTUS ../VIRTUS \
    --genomeDir_human ../VIRTUS/index/STAR_index_human \
    --genomeDir_virus ../VIRTUS/index/STAR_index_virus \
    --figsize 3,3 
```

## Tips

- cwltool may occupy all the system disk by tmp directory. If you suspect the situation, check `/tmp` or avoid by cwltool's option. The example is below. You can also delete the dir by `--rm-tmpdir`.

```
cwltool --tmp-outdir-prefix=/home/yyasumizu/tmp_cwl/ \
--tmpdir-prefix=/home/yyasumizu/tmp_cwl/ \
~/yyoshiaki-git/VIRTUS/workflow/VIRTUS.PE.cwl \
--fastq1 /home/yyasumizu/NGS_public/PRJEB31829_Blimph_EB/donor1_day0_1.fastq.gz \
--fastq2 /home/yyasumizu/NGS_public/PRJEB31829_Blimph_EB/donor1_day0_2.fastq.gz \
--genomeDir_human /home/yyasumizu/yyoshiaki-git/VIRTUS/test/STAR_index_human \
--genomeDir_virus /home/yyasumizu/yyoshiaki-git/VIRTUS/test/STAR_index_virus \
--outFileNamePrefix_human /home/yyasumizu/EB_VIRTUS/donor1_day0/human --nthreads 20
```

- when you specify .cwl files in the absolute path, error may occur. use the relative path.
- note that you cannnot use `\`in --outFileNamePrefix_*
- STAR will require memory at least 30GB. Check your resources.
- Some options instead of docker are available. Specify cwltool option `--user-space-docker-cmd=udocker`, `--singularity`. Increase a limited number by `ulimit -n 4096` (or more) may be required.
- You can specify another host's reference URL such as the mouse in createindex steps, but note that virus references are designed for human viruses. We don't guarantee the result when you changed the reference species.

## VIRTUS2 for single cell RNAseq

### virus detection for 10x or Dropseq

10x and Dropseq use paired-end sequences. The second fastq file contains only transcript's sequences. We recommend users first to run `VIRTUS.SE.cwl` for the second reads　as a virome-wide screening. For detected viruses, users can quantify in several ways; 1) [cellrenger with a modified reference](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr), 2) run [STAR solo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#barcode-and-cdna-on-the-same-mate) with a custom reference, 3) run alevin for the detected virus. `createindex_singlevirus.cwl` in VIRTUS1 can be used for building the index for [Alevin](https://salmon.readthedocs.io/en/latest/alevin.html). For example, the Dropseq's output from SRR8315715 can be screened like the command below.

```
./VIRTUS.SE.cwl \
--fastq ../test/SRR8315715_2.fastq.gz \
--genomeDir_human ../test/STAR_index_human \
--genomeDir_virus ../test/STAR_index_virus \
--outFileNamePrefix_human human \
--nthreads 40
```

### virus detection for SmartSeq2

Just use `VIRTUS.PE.cwl` on each cell individually. When the number of reads is insufficient, VIRTUS may not detect viruses.

## Test

After you clone this repo, you can test VIRTUS2.

```
cd test
bash test.sh
```

## Tool versions

- samtools: 1.15
- bedtools: 2.29.2
- complexity: v0.3.6
- fastp: 0.20.0
- STAR: 2.7.1a

## cwl sources

- [https://github.com/pitagora-network/DAT2-cwl](https://github.com/pitagora-network/DAT2-cwl/tree/develop) : most tools
- [https://github.com/nigyta/bact_genome](https://github.com/nigyta/bact_genome) : fastp

## Contact

Yoshiaki Yasumizu ([yyasumizu@ifrec.osaka-u.ac.jp](yyasumizu@ifrec.osaka-u.ac.jp))
Please use [issue](https://github.com/yyoshiaki/VIRTUS2/issues) first!

## Citation

Yoshiaki Yasumizu, Atsushi Hara, Shimon Sakaguchi, Naganari Ohkura, VIRTUS: a pipeline for comprehensive virus analysis from conventional RNA-seq data, *Bioinformatics*, Volume 37, Issue 10, 15 May 2021, Pages 1465–1467, [https://doi.org/10.1093/bioinformatics/btaa859](https://doi.org/10.1093/bioinformatics/btaa859)

## Acknowledgement

- [Kozo Nishida](https://github.com/kozo2) : CI config and many other supports.
- Ayano Onishi : Logo design.

## Licence

This software is freely available for academic users. Usage for commercial purposes is not allowed. Please refer to the LICENCE page.

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a>
