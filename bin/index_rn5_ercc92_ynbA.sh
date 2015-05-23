#!/usr/bin/env bash

bin/_prepare_goldenPath.sh rn5
bin/_prepare_ercc92.sh
bin/_prepare_ynbA.sh

dir=src/ebwt/rn5_ercc92_ynbA
base=$dir/ref
fa=${base}.fa

mkdir -p $dir
unpigz -c tmp/rn5.fa.gz tmp/ercc92.fa.gz tmp/ynbA.fa.gz > $fa
bowtie-build $fa $base
samtools faidx $fa
