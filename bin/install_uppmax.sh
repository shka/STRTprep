#!/usr/bin/env bash
# [[file:~/Documents/org/Rackham.org::*STRTprep3%20configulation][STRTprep3 configulation:1]]
set -e

. bin/setup_uppmax.sh

gem install bundler
bundle
gem install rake

R --vanilla --quiet <<EOF
source("http://bioconductor.org/biocLite.R")
biocLite(c('devtools', 'impute', 'yaml'), ask=F, lib.loc=.libPaths()[1], lib=.libPaths()[1])
library(devtools)
install_github('shka/samr', ref='test_multblock')
install_github('shka/R-SAMstrt')
EOF

gcc -o bin/_step1b_fastq2oneLine bin/_step1b_fastq2oneLine.c
gcc -o bin/_step1b_fastqs2oneLine bin/_step1b_fastqs2oneLine.c
gcc -o bin/_step1b_trimWithQCFilter bin/_step1b_trimWithQCFilter.c
# STRTprep3 configulation:1 ends here
