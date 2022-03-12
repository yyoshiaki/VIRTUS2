#!/usr/bin/env python3

# %% Package import
import subprocess
import numpy as np
import pandas as pd
import argparse
import os
from scipy import stats
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import matplotlib as mpl
import pathlib
import re
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# %% Argument setting
parser = argparse.ArgumentParser()

parser.add_argument('input_path')
parser.add_argument('--VIRTUSDir', default = os.path.dirname(os.path.abspath(__file__)))
parser.add_argument('--genomeDir_human', required = True)
parser.add_argument('--genomeDir_virus', required = True)
parser.add_argument('--outFileNamePrefix_human', default = 'human')
parser.add_argument('--nthreads', default = '16')
parser.add_argument('-s', '--Suffix_SE')
parser.add_argument('-s1', '--Suffix_PE_1')
parser.add_argument('-s2', '--Suffix_PE_2')
parser.add_argument('--fastq', action = 'store_true')

args = parser.parse_args()

# %%
df = pd.read_csv(args.input_path)
df.columns = ["Name", "SRR", "Layout", "Group"] + list(df.columns[4:])
first_dir = os.getcwd()

print(args.VIRTUSDir)
try:
    if os.path.exists(os.path.join(args.VIRTUSDir, "workflow/VIRTUS.PE.cwl")):
        dir_VIRTUS = os.path.join(args.VIRTUSDir, "workflow")
    elif os.path.exists(os.path.join(args.VIRTUSDir, "VIRTUS.PE.cwl")):  
        dir_VIRTUS = args.VIRTUSDir
    else:
        raise ValueError('not found VIRTUS.PE.cwl or VIRTUS.SE.cwl')
except (ValueError, IndexError):
    exit('invalid path to VIRTUS. try to change --VIRTUSDir to the absolute path.')

# %%
series_list = []

for index, item in df.iterrows():
    # parameter setting
    if args.fastq == True:
        name = item['Name']
        # dir = os.path.dirname(item["SRR"])
        sample_index = os.path.basename(item["SRR"])
        # os.chdir(dir)
        p_temp = pathlib.Path(".")
        files = [str(i) for i in list(p_temp.iterdir()) if i.is_file()]
        if item["Layout"] == "PE":
            if args.Suffix_PE_1 == None:
                pattern_1 = "^" + sample_index + "_1((\.fq\.gz)|(\.fq)|(\.fastq)|(\.fastq\.gz))$"
                matched_files_1 = sorted([i for i in files if re.match(pattern_1,i)])
                if len(matched_files_1) >= 1:
                    fastq1 = matched_files_1[0]
                    print("fastq_1:",fastq1)
                else:
                    print("fastq_1 not found")
            else:
                fastq1 = sample_index + args.Suffix_PE_1
            if args.Suffix_PE_2 == None:
                pattern_2 = "^" + sample_index + "_2((\.fq\.gz)|(\.fq)|(\.fastq)|(\.fastq\.gz))$"
                matched_files_2 = sorted([i for i in files if re.match(pattern_2,i)])
                if len(matched_files_2) >= 1:
                    fastq2 = matched_files_2[0]
                    print("fastq_2:",fastq2)
                else:
                    print("fastq_2 not found")
            else:
                fastq2 = sample_index + args.Suffix_PE_2
        elif item["Layout"] == "SE":
            if args.Suffix_SE == None:
                pattern = "^" + sample_index + "((\.fq\.gz)|(\.fq)|(\.fastq)|(\.fastq\.gz))$"
                matched_files = sorted([i for i in files if re.match(pattern,i)])
                if len(matched_files) >= 1:
                    fastq = matched_files[0]
                    print("fastq:",fastq)
                else:
                    print("fastq not found")
            else:
                fastq = sample_index + args.Suffix_SE
        else:
            print("Layout Error")
    else:
        name = item['Name']
        dir = item["SRR"]
        sample_index = item["SRR"]
        prefetch_cmd = " ".join(["prefetch",sample_index])
        fasterq_cmd = " ".join(["fasterq-dump", "--split-files", sample_index + ".sra", "-e","16"])
        if item["Layout"] == "PE":
            fastq1 = sample_index + "_1.fastq.gz"
            fastq2 = sample_index + "_2.fastq.gz"
            input_list = [fastq1,fastq2]
        elif item["Layout"] == "SE":
            fastq = sample_index + ".fastq.gz"
            input_list = [fastq]
        else:
            print("Layout Error")

    if item["Layout"] =="PE":
        VIRTUS_cmd = " ".join([
            "cwltool",
            "--rm-tmpdir",
            os.path.join(dir_VIRTUS, "VIRTUS.PE.cwl"), 
            "--fastq1", '../'+fastq1,
            "--fastq2", '../'+fastq2, 
            "--genomeDir_human", args.genomeDir_human, 
            "--genomeDir_virus", args.genomeDir_virus,
            "--outFileNamePrefix_human", args.outFileNamePrefix_human,
            "--nthreads", args.nthreads,
            "--filename_output", "VIRTUS.{}.txt".format(name)
        ])
    elif item["Layout"] =="SE":
        VIRTUS_cmd = " ".join([
            "cwltool",
            "--rm-tmpdir",
            os.path.join(dir_VIRTUS, "VIRTUS.SE.cwl"), 
            "--fastq", '../'+fastq,
            "--genomeDir_human", args.genomeDir_human, 
            "--genomeDir_virus", args.genomeDir_virus,
            "--outFileNamePrefix_human", args.outFileNamePrefix_human,
            "--nthreads", args.nthreads,
            "--filename_output", "VIRTUS.{}.txt".format(name)
        ])
    else:
        print("Layout Error")

    # run
    os.makedirs(name, exist_ok=True)

    if args.fastq == False:
        if item["Layout"] =="PE":
            is_exit_fq = os.path.exists(dir+'/'+fastq1)
        else:
            is_exit_fq = os.path.exists(dir+'/'+fastq)
        if not is_exit_fq:
            print(prefetch_cmd,"\n")    
            subprocess.run(prefetch_cmd, shell = True)

            print(fasterq_cmd,"\n")
            subprocess.run(fasterq_cmd, shell = True)
            if item["Layout"] == "PE":
                os.rename(sample_index + ".sra_1.fastq.gz", fastq1)
                os.rename(sample_index + ".sra_2.fastq.gz", fastq2)
            elif item["Layout"] == "SE":
                os.rename(sample_index + ".sra.fastq", fastq)
            
            pigz_cmd = "pigz *.fastq "
            print(pigz_cmd, "\n")
            subprocess.run(pigz_cmd, shell = True)

    print(VIRTUS_cmd,"\n")
    os.chdir(name)
    subprocess.run(VIRTUS_cmd, shell = True)
        
    if args.fastq == False:
        for i in input_list:
            pigz_cmd = " ".join(["pigz", i])
            print(pigz_cmd, "\n")
            subprocess.run(pigz_cmd, shell = True)
            
    df_virus = pd.read_table("VIRTUS.{}.txt".format(name), index_col = 0)
    series_virus = df_virus.loc[:,"rate_hit"]
    series_virus = series_virus.rename(item['Name'])
    series_list.append(series_virus)

    os.chdir("..")

# %% summary
summary = pd.concat(series_list, axis = 1).fillna(0).T
summary["Group"] = df["Group"].values

summary_dict = {}
Group = summary["Group"].unique()
for i in Group:
    summary_dict[i] = summary[summary["Group"] == i]

uval = pd.Series(dtype = "float64")
pval = pd.Series(dtype = "float64")

if summary["Group"].nunique() == 2:
    print("Conducting Mann-Whitney U-test")
    for i in range(0,len(summary.columns)-1):
        if summary["Group"].nunique() == 2:

            u, p = stats.mannwhitneyu(summary_dict[Group[0]].iloc[:,i],summary_dict[Group[1]].iloc[:,i], alternative = "two-sided")
            uval[summary.columns[i]] = u
            pval[summary.columns[i]] = p

    fdr = pd.Series(multipletests(pval,method = "fdr_bh")[1], index = pval.index)

    summary.loc["u-value"] = uval
    summary.loc["p-value"] = pval
    summary.loc["FDR"] = fdr

summary.to_csv("summary.csv")

# %% Graph drawing
g = sns.clustermap(summary.iloc[:-3,:-1].T, method = "ward", metric="euclidean")
g.savefig("clustermap.pdf", bbox_inches='tight')

print('All processes succeeded.')
