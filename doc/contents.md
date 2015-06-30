# STRTprep

STRTprep (https://github.com/shka/STRTprep) is an open-source software package for preprocess and analysis of STRT RNA-seq data. The first version was for the preprocessing and quality-check purposes, and the version 2 contains automated statistical tests of the differential expression. This package can run on both Linux and OSX.

Edit two files to describe your experiment and study design, and run the pipeline as below, then you can finish from preprocess of raw reads until differential expression tests. Enjoy more further downstream analysis!

```bash
# STEP 1: Installation
## pause; decide name and the location of your project - now './STRTprep2.test'
git clone https://github.com/shka/STRTprep.git STRTprep2.test
cd STRTprep2.test
bin/install.sh
. bin/setup.sh
## pause; choose reference genome and transcriptome - now hg19 and RefSeq
bin/index_hg19_ercc92_ynbA_u13369.sh
bin/index_refGene.sh hg19 src/ebwt/hg19_ercc_ynbA_u13369/ref

# STEP 2: Design of your experiments
## pause; edit conf.yaml and src/samples.csv

# STEP 3: Preprocessing & quality check
rake -m -j `gnproc` qc > qc.log 2>&1
## pause; check qc.log & out/byGene/samples.xls

# STEP 4: Differential expression analysis
rake gene > gene.log 2>&1
## pause; check gene.log, out/byGene/diffexp.xls & out/byGene/plugins_*
##        go back to the step 2 when you would like to compare more

# DONE! Enjoy more further analysis!
```

> This document is for the version 2, but version 3 beta (branch `v3dev`) is ongoing; the version 3 will provide TFE-based quantification to analyze novel genes/start-sites including long-noncoding RNAs.

## Table of contents

1. [Installation of STRTprep](install.md)
2. [Protocol of preprocessing and analysis](protocol.md)
3. [Interpretation of the results](result.md)

## Notation

> Sentences with a vertical box at the left like this are usually tips or pitfalls.

```bash
Sentences within a box are codes; usually these are commands to be executed on your bash console.
```
