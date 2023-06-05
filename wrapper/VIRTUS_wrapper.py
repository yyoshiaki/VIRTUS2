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
import matplotlib.pyplot as plt
import pathlib
import re
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

from scattermap import scattermap

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
parser.add_argument('--figsize', default = '8,3', help='(default:8,3)')
parser.add_argument('--th_cov', default = '10', help='threshold of max viral coverage to plot, test (default:10)')
parser.add_argument('--th_rate', default = '0.0001', help='threshold of max rate virus/human to plot, test (default:0.0001)')

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
list_df_res = []

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
        fasterq_cmd = " ".join(["fasterq-dump", "--split-files", sample_index + '/' + sample_index + ".sra", "-e","16"])
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
            "--filename_output", "VIRTUS.output.txt"
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
            "--filename_output", "VIRTUS.output.txt"
        ])
    else:
        print("Layout Error")

    # run
    os.makedirs(name, exist_ok=True)

    if args.fastq == False:
        if item["Layout"] =="PE":
            is_exit_fq = os.path.exists(fastq1)
        else:
            is_exit_fq = os.path.exists(fastq)
        if not is_exit_fq:
            print(prefetch_cmd,"\n")    
            subprocess.run(prefetch_cmd, shell = True)

            print(fasterq_cmd,"\n")
            subprocess.run(fasterq_cmd, shell = True)

            pigz_cmd = "pigz *.fastq "
            print(pigz_cmd, "\n")
            subprocess.run(pigz_cmd, shell = True)

            # if item["Layout"] == "PE":
            #     os.rename(sample_index + ".sra_1.fastq.gz", fastq1)
            #     os.rename(sample_index + ".sra_2.fastq.gz", fastq2)
            # elif item["Layout"] == "SE":
            #     os.rename(sample_index + ".sra.fastq", fastq)
            

    print(VIRTUS_cmd,"\n")
    os.chdir(name)
    if not os.path.exists('VIRTUS.output.txt'):
        subprocess.run(VIRTUS_cmd, shell = True)
        
            
    d = pd.read_table("VIRTUS.output.txt".format(name))
    list_df_res.append(d)
    d['name'] = item['Name']

    os.chdir("..")

# %% summary
df_res = pd.concat(list_df_res)

th_cov = float(args.th_cov)
th_rate = float(args.th_rate)

df_rate = pd.pivot(df_res, index='virus', columns='name', values="rate_hit").fillna(0)
df_cov = pd.pivot(df_res, index='virus', columns='name', values="coverage").fillna(0)

df_cov = df_cov[df_cov.max(axis=1) > th_cov]
df_rate = df_rate[df_rate.max(axis=1) > th_rate]

list_index = (set(list(df_rate.index)) & set(list(df_cov.index)))

df_cov = df_cov.reindex(list_index)
df_rate = df_rate.reindex(list_index)

df_cov = df_cov[df_rate.columns]

summary = pd.merge(df_cov, df_rate, left_index=True, right_index=True, suffixes=('_cov', '_rate'))


uval = pd.Series(dtype = "float64")
pval = pd.Series(dtype = "float64")

# print(summary)

print("### Stats ###")
print('Threshold coverage : ', th_cov)
print('Threshold rate : ', th_rate)

print('# detected viruses : ', df_res.shape[0])
print('Max coverage : ', df_res.coverage.max())
print('Max rate : ', df_res.rate_hit.max())

print('# retained viruses :', summary.shape[0])

if df["Group"].nunique() == 2:
    if summary.shape[0] > 0:
        print("Conducting Mann-Whitney U-test")
        for v in summary.index:
            u, p = stats.mannwhitneyu(summary.loc[v, ['{}_rate'.format(x) for x in df.loc[df['Group']==df['Group'].unique()[0], 'Name']]],
            summary.loc[v, ['{}_rate'.format(x) for x in df.loc[df['Group']==df['Group'].unique()[1], 'Name']]], alternative = "two-sided")
            uval[v] = u
            pval[v] = p

        fdr = pd.Series(multipletests(pval,method = "fdr_bh")[1], index = pval.index)
    else:
        print("Skipped Mann-Whitney U-test")
        fdr = ""

    summary["u-value"] = uval
    summary["p-value"] = pval
    summary["FDR"] = fdr

summary.to_csv("summary.csv")


# %% Graph drawing
if summary.shape[0] > 0:
    figsize = (int(args.figsize.split(',')[0]), int(args.figsize.split(',')[1]))
    # g = sns.clustermap(summary, method = "ward", metric="euclidean", figsize=figsize)
    # g.savefig("clustermap.pdf", bbox_inches='tight')


    with sns.axes_style("white"):
        plt.figure(figsize=figsize)
        ax = scattermap(df_rate, square=True, marker_size=df_cov, cmap='viridis_r',
        cbar_kws={'label': 'v/h rate'})

        #make a legend:
        pws = [20, 40, 60, 80, 100]
        for pw in pws:
            plt.scatter([], [], s=(pw), c="k",label=str(pw))

        h, l = plt.gca().get_legend_handles_labels()
        plt.legend(h[1:], l[1:], labelspacing=.3, title="coverage(%)", borderpad=0, framealpha=0, edgecolor="w",
                    bbox_to_anchor=(1.1, -.1), ncol=1, loc='upper left', borderaxespad=0)
        plt.savefig('scattermap.pdf' , bbox_inches='tight')

print('All processes succeeded.')
