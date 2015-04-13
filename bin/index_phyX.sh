#!/usr/bin/env bash

seq=src/NC_001422.fa

while [ -z "`grep '^>' $seq`" ]; do
    curl -o $seq "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=9626372&strand=1&rettype=fasta&retmode=text"
done

dir=src/ebwt/phyX
base=$dir/ref
ref=${base}.fa

mkdir -p $dir
gawk '/^>/{ print ">NC_001422.1" } !/^>/{ printf $1 }' $seq | gfold -w 50 > $ref
bowtie-build $ref $base
samtools faidx $ref
