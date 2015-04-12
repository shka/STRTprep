#!/usr/bin/env bash

curl -o src/canFam3.fa.gz 'http://hgdownload.cse.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.fa.gz'

bin/_prepare_ercc92.sh
bin/_prepare_ynbA.sh

dir=src/ebwt/canFam3_ercc92_ynbA
base=$dir/ref
fa=${base}.fa

mkdir -p $dir
unpigz -c src/canFam3.fa.gz tmp/ercc92.fa.gz tmp/ynbA.fa.gz > $fa
bowtie-build $fa $base
samtools faidx $fa
