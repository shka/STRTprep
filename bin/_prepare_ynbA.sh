#!/usr/bin/env bash

seq=src/EF011072.fa

while [ -z "`grep '^>' $seq`" ]; do
    curl -o $seq "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=116733912&strand=1&rettype=fasta&retmode=text"
done

mkdir -p tmp
gawk '/^>/{ print ">RNA_SPIKE_ynbA" } /^ACC/{ sub(/^AC/, "GGGGGAATT"); printf $1 } !/^(>|ACC)/{ printf $1 } END{ print }' $seq | gfold -w 50 | gzip -c > tmp/ynbA.fa.gz
