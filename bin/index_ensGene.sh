#!/usr/bin/env bash

ver=$1
dir=src/ebwt/${ver}_ensGene
base=$dir/ref
gtf=${base}.gtf

mkdir -p $dir
curl -o $dir/ensGene.txt.gz "http://hgdownload.cse.ucsc.edu/goldenPath/$ver/database/ensGene.txt.gz"
curl -o $dir/ensemblToGeneName.txt.gz "http://hgdownload.cse.ucsc.edu/goldenPath/$ver/database/ensemblToGeneName.txt.gz"

gunzip -c $dir/ensGene.txt.gz | gcut -f 2-11 | genePredToGtf file stdin $gtf
bin/_preprocess_annotation_ensGene.rb $dir/ensemblToGeneName.txt.gz $dir/ensGene.txt.gz $base

tophat -p `gnproc` --bowtie1 -G $gtf --transcriptome-index $base $2
