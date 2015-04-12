#!/usr/bin/env bash

curl -o src/U13369.fa "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=555853&strand=1&rettype=fasta&retmode=text"

mkdir -p tmp
gawk '/^>/{ print ">RIBO_U13369.1" } !/^>/{ printf $1 }' src/U13369.fa | gfold -w 50 | gzip -c > tmp/u13369.fa.gz
