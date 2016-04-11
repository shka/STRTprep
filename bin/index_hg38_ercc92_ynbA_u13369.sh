#!/usr/bin/env bash

set -e

bin/_prepare_goldenPath.sh hg38
bin/_prepare_ercc92.sh
bin/_prepare_ynbA.sh
bin/_prepare_u13369.sh

if [ -z $1 ]; then
    dir=src/ebwt/hg38_ercc92_ynbA_u13369
else
    dir=$1
fi

base=$dir/ref
fa=${base}.fa

mkdir -p $dir
unpigz -c tmp/hg38.fa.gz tmp/ercc92.fa.gz tmp/ynbA.fa.gz tmp/u13369.fa.gz > $fa
if [ -n $2 ]; then
    cat $2 | gfold -w 50 >> $fa
fi
bowtie-build $fa $base
samtools faidx $fa
