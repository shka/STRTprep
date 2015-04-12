#!/usr/bin/env bash

curl -o src/NC_001422.fa "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=9626372&strand=1&rettype=fasta&retmode=text"

mkdir -p src/ebwt/phyX
gawk '/^>/{ print ">NC_001422.1" } !/^>/{ printf $1 }' src/NC_001422.fa | gfold -w 50 > src/ebwt/phyX/ref.fa
bowtie-build src/ebwt/phyX/ref.fa src/ebwt/phyX/ref
samtools faidx src/ebwt/phyX/ref.fa
