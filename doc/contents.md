# STRTprep

STRTprep (https://github.com/shka/STRTprep) is an open-source software package for preprocess and analysis of STRT RNA-seq data. The first version was for the preprocessing and quality-check purposes. The version 2 realized automated statistical tests of the differential expression. The version 3 performs "transcription start region" based high resolution analysis including long noncoding RNAs and novel genes/transcripts. Then finally, the version 4 is an optimization of complex studies. This package can run on both Linux and OSX.

Edit `conf.yaml` to describe your experiment and study design, and run the pipeline as below, then you can finish from preprocess of raw reads until differential expression tests. Enjoy more further downstream analysis!

```bash
git clone -b v4dev https://github.com/shka/STRTprep.git STRTprep4.test
cd STRTprep4.test
bin/install.sh
. bin/setup.sh
```

> This document is for the version 4 beta (branch `v4dev`).

## Notation

> Sentences with a vertical box at the left like this are usually tips or pitfalls.

```bash
Sentences within a box are codes; usually these are commands to be executed on your bash console.
```
