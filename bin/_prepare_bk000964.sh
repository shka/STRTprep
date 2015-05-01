#!/usr/bin/env bash

seq=src/BK000964.fa

while [ -z "`grep '^>' $seq`" ]; do
    curl -o $seq "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=511668571&strand=1&rettype=fasta&retmode=text"
done

mkdir -p tmp
gawk '/^>/{ print ">RIBO_BK000964" } !/^>/{ printf $1 } END{ print }' $seq | gfold -w 50 | gzip -c > tmp/bk000964.fa.gz
