#!/usr/bin/env bash

set -e

. bin/setup_uppmax.sh

if ! [ -d '.homebrew' ]; then
    git clone https://github.com/Homebrew/linuxbrew.git .homebrew
fi
## Issue: Unable to bootstrap gcc
##   https://github.com/Homebrew/linuxbrew/issues/137
ln -s `which gcc` \
   `brew --prefix`/bin/gcc-`gcc -dumpversion | cut -d. -f1,2` && true
ln -s `which g++` \
   `brew --prefix`/bin/g++-`g++ -dumpversion | cut -d. -f1,2` && true
ln -s `which gfortran` \
   `brew --prefix`/bin/gfortran-`gfortran -dumpversion | cut -d. -f1,2` && true

brew tap homebrew/science
brew prune
brew update

gem install bundler
bundle
gem install rake

R --vanilla --quiet <<EOF
source("http://bioconductor.org/biocLite.R")
biocLite(c('devtools', 'samr', 'yaml'), ask=F, lib.loc=.libPaths()[1], lib=.libPaths()[1])
library(devtools)
install_github('shka/samr', ref='test_multblock')
install_github('shka/R-SAMstrt')
EOF

gcc -o bin/_step1b_fastq2oneLine bin/_step1b_fastq2oneLine.c
gcc -o bin/_step1b_fastqs2oneLine bin/_step1b_fastqs2oneLine.c
gcc -o bin/_step1b_trimWithQCFilter bin/_step1b_trimWithQCFilter.c

