#!/usr/bin/env bash

set -e
. bin/setup.sh

sudo apt-get install build-essential curl file git zlib1g-dev emacs-nox

sh -c "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh)"
test -d ~/.linuxbrew && eval $(~/.linuxbrew/bin/brew shellenv)
test -d /home/linuxbrew/.linuxbrew && eval $(/home/linuxbrew/.linuxbrew/bin/brew shellenv)
test -r ~/.bash_profile && echo "eval \$($(brew --prefix)/bin/brew shellenv)" >>~/.bash_profile
echo "eval \$($(brew --prefix)/bin/brew shellenv)" >>~/.profile

brew tap brewsci/bio
brew tap brewsci/science
brew install bedtools bowtie coreutils kent-tools parallel pigz r ruby samtools@0.1

wget -O src/bowtie-1.2.3-linux-x86_64.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.3/bowtie-1.2.3-linux-x86_64.zip/download
unzip src/bowtie-1.2.3-linux-x86_64.zip -d src/
ln -s $PWD/src/bowtie-1.2.3-linux-x86_64/bowtie* bin/

wget -O src/tophat-2.1.1.Linux_x86_64.tar.gz https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
tar zxvf src/tophat-2.1.1.Linux_x86_64.tar.gz -C src/
ln -s `find $PWD/src/tophat-2.1.1.Linux_x86_64/ -type f -maxdepth 1 -executable` bin/

wget -O src/cufflinks-2.2.1.Linux_x86_64.tar.gz http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar zxvf src/cufflinks-2.2.1.Linux_x86_64.tar.gz -C src/
ln -s `find $PWD/src/cufflinks-2.2.1.Linux_x86_64/ -type f -maxdepth 1 -executable` bin/

bundle

R --vanilla --quiet <<EOF
## source("http://bioconductor.org/biocLite.R")
install.packages("BiocManager")
BiocManager::install(c('remotes', 'yaml', 'impute', 'RColorBrewer', 'beeswarm'))
library(remotes)
install_github('renozao/pkgmaker@develop')
install_github('renozao/NMF@devel')
install_github('shka/samr', ref='test_multblock')
install_github('shka/R-SAMstrt', dependencies=F, upgrade='never')
EOF

gcc -o bin/_step1b_fastq2oneLine bin/_step1b_fastq2oneLine.c
gcc -o bin/_step1b_fastqs2oneLine bin/_step1b_fastqs2oneLine.c
gcc -o bin/_step1b_trimWithQCFilter bin/_step1b_trimWithQCFilter.c
