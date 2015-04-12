#!/usr/bin/env bash

bin/_prepare_goldenPath.sh hg19
bin/_prepare_ercc92.sh
bin/_prepare_ynbA.sh
bin/_prepare_u13369.sh

dir=src/ebwt/hg19_ercc92_ynbA_u13369
base=$dir/ref
fa=${base}.fa

mkdir -p $dir
unpigz -c tmp/hg19.fa.gz tmp/ercc92.fa.gz tmp/ynbA.fa.gz tmp/u13369.fa.gz > $fa
bowtie-build $fa $base
samtools faidx $fa
