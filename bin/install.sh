#!/usr/bin/env bash

set -e
. bin/setup.sh

sudo apt-get install build-essential curl file git zlib1g-dev

sh -c "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh)"
test -d ~/.linuxbrew && eval $(~/.linuxbrew/bin/brew shellenv)
test -d /home/linuxbrew/.linuxbrew && eval $(/home/linuxbrew/.linuxbrew/bin/brew shellenv)
test -r ~/.bash_profile && echo "eval \$($(brew --prefix)/bin/brew shellenv)" >>~/.bash_profile
echo "eval \$($(brew --prefix)/bin/brew shellenv)" >>~/.profile

brew tap brewsci/bio
brew tap brewsci/science
brew install bedtools coreutils kent-tools parallel pigz r ruby samtools@0.1

mkdir -p ~/.local/bin

wget -O src/bowtie-1.1.2-linux-x86_64.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip/download
unzip src/bowtie-1.1.2-linux-x86_64.zip -d ~/.local/
ln -s ~/.local/bowtie-1.1.2-linux-x86_64/bowtie* ~/.local/bin/

wget -O src/tophat-2.1.1.Linux_x86_64.tar.gz https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
tar zxvf src/tophat-2.1.1.Linux_x86_64.tar.gz -C ~/.local/
ln -s `find ~/.local/tophat-2.1.1.Linux_x86_64/ -maxdepth 1 -type f -executable` ~/.local/bin/

wget -O src/cufflinks-2.2.1.Linux_x86_64.tar.gz http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar zxvf src/cufflinks-2.2.1.Linux_x86_64.tar.gz -C ~/.local/
ln -s `find ~/.local/cufflinks-2.2.1.Linux_x86_64/ -maxdepth 1 -type f -executable` ~/.local/bin/

bundle

R --vanilla --quiet <<EOF
install.packages("BiocManager", repos="http://cran.r-project.org")
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
