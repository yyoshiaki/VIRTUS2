#!/bin/bash
set -xe

DIR_VIRTUS=$HOME/Programs/VIRTUS2
DIR_INDEX_ROOT=$HOME/reference/VIRTUS_2.0

python3 VIRTUS_wrapper.py input.fastq.csv \
    --VIRTUSDir $DIR_VIRTUS \
    --genomeDir_human $DIR_INDEX_ROOT/STAR_index_human \
    --genomeDir_virus $DIR_INDEX_ROOT/STAR_index_virus \
    --nthreads=4 --fastq
