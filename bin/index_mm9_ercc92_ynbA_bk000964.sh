#!/usr/bin/env bash

bin/_prepare_goldenPath.sh mm9
bin/_prepare_ercc92.sh
bin/_prepare_ynbA.sh
bin/_prepare_bk000964.sh

dir=src/ebwt/mm9_ercc92_ynbA_bk000964
base=$dir/ref
fa=${base}.fa

mkdir -p $dir
unpigz -c tmp/mm9.fa.gz tmp/ercc92.fa.gz tmp/ynbA.fa.gz tmp/bk000964.fa.gz > $fa
bowtie-build $fa $base
samtools faidx $fa
