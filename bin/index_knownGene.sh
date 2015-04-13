#!/usr/bin/env bash

ver=$1
dir=src/ebwt/${ver}_knownGene
base=$dir/ref
gtf=${base}.gtf

mkdir -p $dir
curl -o $dir/knownGene.txt.gz "http://hgdownload.cse.ucsc.edu/goldenPath/$ver/database/knownGene.txt.gz"
curl -o $dir/kgXref.txt.gz "http://hgdownload.cse.ucsc.edu/goldenPath/$ver/database/kgXref.txt.gz"

gunzip -c $dir/knownGene.txt.gz | gcut -f 1-10 | genePredToGtf file stdin $gtf
bin/_preprocess_annotation_knownGene.rb $dir/kgXref.txt.gz $dir/knownGene.txt.gz $base

tophat -p `gnproc` --bowtie1 -G $gtf --transcriptome-index $base $2
