#!/usr/bin/env bash

curl -o src/danRer7.fa.gz 'http://hgdownload.cse.ucsc.edu/goldenPath/danRer7/bigZips/danRer7.fa.gz'

bin/_prepare_ercc92.sh
bin/_prepare_ynbA.sh

dir=src/ebwt/danRer7_ercc92_ynbA
base=$dir/ref
fa=${base}.fa

mkdir -p $dir
unpigz -c src/danRer7.fa.gz tmp/ercc92.fa.gz tmp/ynbA.fa.gz > $fa
bowtie-build $fa $base
samtools faidx $fa
