#!/usr/bin/env bash

curl -o src/ERCC92.zip 'https://tools.lifetechnologies.com/content/sfs/manuals/ERCC92.zip'

mkdir -p tmp
unzip -cq src/ERCC92.zip ERCC92.fa | gawk 'BEGIN{p=""} /^>/{ print p ">RNA_SPIKE_" substr($1, 2); printf("GGGGGAATTC") } /^[ACGT]/{ printf $1; p="\n" } END{ print }' | gfold -w 50 | gzip -c > tmp/ercc92.fa.gz
