#!/usr/bin/env bash

ver=$1
dir=src/ebwt/${ver}_refGene
base=$dir/ref
gtf=${base}.gtf

mkdir -p $dir
curl -o $dir/refGene.txt.gz "http://hgdownload.cse.ucsc.edu/goldenPath/$ver/database/refGene.txt.gz"

gunzip -c $dir/refGene.txt.gz | gcut -f 2-11 | genePredToGtf file stdin $gtf
bin/_preprocess_annotation_refGene.rb $dir/refGene.txt.gz $dir/refGene.txt.gz $base

tophat -p `gnproc` --bowtie1 -G $gtf --transcriptome-index $base $2
