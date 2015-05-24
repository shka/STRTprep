#!/usr/bin/env bash

## Issues:
##

set -e
. bin/setup.sh

if ! [ -d '.homebrew' ]; then
    if [ `uname` = 'Linux' ]; then
	git clone https://github.com/Homebrew/linuxbrew.git .homebrew
    else
	mkdir .homebrew
	curl -L https://github.com/Homebrew/homebrew/tarball/master | tar xz --strip 1 -C .homebrew
    fi
    ## Issue: Unable to bootstrap gcc
    ##   https://github.com/Homebrew/linuxbrew/issues/137
    ln -s `which gcc` \
       `brew --prefix`/bin/gcc-`gcc -dumpversion | cut -d. -f1,2` && true
    ln -s `which g++` \
       `brew --prefix`/bin/g++-`g++ -dumpversion | cut -d. -f1,2` && true
    ln -s `which gfortran` \
       `brew --prefix`/bin/gfortran-`gfortran -dumpversion | cut -d. -f1,2` && true
fi

brew update
brew tap homebrew/science
brew install coreutils

## Issue: Bedtools not later than 2.22.0
##   https://github.com/arq5x/bedtools2/issues/212
if ! [ -d 'bin/bedtools-2.22.0' ]; then
    base=bin/bedtools-2.22.0
    mkdir -p $base
    curl -s https://codeload.github.com/arq5x/bedtools2/tar.gz/v2.22.0 | tar -C $base --strip-component 1 -zxvf -
    cd $base
    make -j `gnproc`
    cd ../..
fi

## Issue: Hidden dependency in R 3.2.0
brew install curl --build-from-source
brew install R --with-openblas --without-tcltk --without-x11

brew install ruby
brew install parallel
brew install pigz
brew install gawk
brew install samtools-0.1
brew install bowtie
brew install tophat --with-bowtie
brew install kent-tools

## Issue: Cufflinks 2.2.1 is unstable
##   https://github.com/Homebrew/homebrew-science/issues/2254
if ! [ -d 'bin/cufflinks-2.1.1' ]; then
    mkdir -p bin/cufflinks-2.1.1
    if [ `uname` = 'Linux' ]; then
	curl -s http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.1.1.Linux_x86_64.tar.gz | tar -C bin/cufflinks-2.1.1 --strip-component 1 -zxvf -
    else
	curl -s http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.1.1.OSX_x86_64.tar.gz | tar -C bin/cufflinks-2.1.1 --strip-component 1 -zxvf -
    fi
fi

gem install bundler
bundle

R --vanilla --quiet <<EOF
source("http://bioconductor.org/biocLite.R")
biocLite(c('devtools', 'samr', 'yaml'), ask=F)
library(devtools)
install_github('shka/samr', ref='test_multblock')
install_github('shka/R-SAMstrt')
EOF

gcc -o bin/_step1b_fastq2oneLine bin/_step1b_fastq2oneLine.c
gcc -o bin/_step1b_trimWithQCFilter bin/_step1b_trimWithQCFilter.c
