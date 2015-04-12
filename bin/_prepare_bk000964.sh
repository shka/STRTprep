#!/usr/bin/env bash

curl -o src/BK000964.fa "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=511668571&strand=1&rettype=fasta&retmode=text"

mkdir -p tmp
gawk '/^>/{ print ">RIBO_BK000964" } !/^>/{ printf $1 }' src/BK000964.fa | gfold -w 50 | gzip -c > tmp/bk000964.fa.gz
