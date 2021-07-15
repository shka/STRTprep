#!/usr/bin/env bash

bin/_prepare_goldenPath.sh hg19
bin/_prepare_ercc92.sh
bin/_prepare_ynbA.sh
bin/_prepare_u13369.sh

if [ -z $1 ]; then
    dir=src/ebwt/hg19_ercc92_ynbA_u13369
else
    dir=$1
fi
base=$dir/ref
fa=${base}.fa

mkdir -p $dir
unpigz -c tmp/hg19.fa.gz tmp/ercc92.fa.gz tmp/ynbA.fa.gz tmp/u13369.fa.gz > $fa
if ! [ -z ${2+UNDEF} ]; then
    cat $2 | gfold -w 50 >> $fa
fi
bowtie-build $fa $base
samtools faidx $fa
