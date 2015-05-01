#!/usr/bin/env bash

curl -o src/susScr3.fa.gz 'http://hgdownload.cse.ucsc.edu/goldenPath/susScr3/bigZips/susScr3.fa.gz'

bin/_prepare_ercc92.sh
bin/_prepare_ynbA.sh

dir=src/ebwt/susScr3_ercc92_ynbA
base=$dir/ref
fa=${base}.fa

mkdir -p $dir
unpigz -c src/susScr3.fa.gz tmp/ercc92.fa.gz tmp/ynbA.fa.gz > $fa
bowtie-build $fa $base
samtools faidx $fa
