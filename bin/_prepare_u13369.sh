#!/usr/bin/env bash

seq=src/U13369.fa

while [ -z "`grep '^>' $seq`" ]; do
    curl -o $seq "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=555853&strand=1&rettype=fasta&retmode=text"
done

mkdir -p tmp
gawk '/^>/{ print ">RIBO_U13369.1" } !/^>/{ printf $1 } END{ print }' $seq | gfold -w 50 | gzip -c > tmp/u13369.fa.gz
