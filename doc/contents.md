# STRTprep

STRTprep (https://github.com/shka/STRTprep) is an open-source software package for preprocess and analysis of STRT RNA-seq data. The first version was for the preprocessing and quality-check purposes, then the version 2 realized automated statistical tests of the differential expression, and finally the version 3 performs "transcription start region" based high resolution analysis including long noncoding RNAs and novel genes/transcripts. This package can run on both Linux and OSX.

Edit two files to describe your experiment and study design, and run the pipeline as below, then you can finish from preprocess of raw reads until differential expression tests. Enjoy more further downstream analysis!

```bash
git clone https://github.com/shka/STRTprep.git STRTprep3.test
cd STRTprep3.test
bin/install.sh
. bin/setup.sh
bin/index_hg19_ercc92_ynbA_u13369.sh
bin/index_refGene.sh hg19 src/ebwt/hg19_ercc_ynbA_u13369/ref
# edit src/conf.yaml & src/samples.csv
rake -m -j `gnproc` qc > qc.log 2>&1
rake gene > gene.log 2>&1
rake > tfe.log 2>&1
```

> This document is for the version 3 beta (branch `v3dev`).

## Table of contents

1. [Installation of STRTprep](install.md)
2. [Protocol of preprocessing and analysis](protocol.md)
3. [Interpretation of the results](result.md)
4. STRTprepHelper - R interface to access the results for the further analysis
5. [Plugin framework](plugin.md) and extension of STRTprep

## Notation

> Sentences with a vertical box at the left like this are usually tips or pitfalls.

```bash
Sentences within a box are codes; usually these are commands to be executed on your bash console.
```
