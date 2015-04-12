#!/usr/bin/env bash

curl -o src/EF011072.fa "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=116733912&strand=1&rettype=fasta&retmode=text"

mkdir -p tmp
gawk '/^>/{ print ">RNA_SPIKE_ynbA" } /^ACC/{ sub(/^AC/, "GGGGGAATT"); printf $1 } !/^(>|ACC)/{ printf $1 }' src/EF011072.fa | gfold -w 50 | gzip -c > tmp/ynbA.fa.gz
