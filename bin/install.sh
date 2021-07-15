#!/usr/bin/env bash

# sudo apt-get install build-essential curl file git procps wget zlib1g-dev

set -e

if [ ! -d "$PWD/homebrew" ]; then
    mkdir $PWD/homebrew \
	&& curl -L https://github.com/Homebrew/brew/tarball/master \
	| tar xz --strip 1 -C $PWD/homebrew
fi

. bin/setup.sh

brew install coreutils
HOMEBREW_MAKE_JOBS=`gnproc` brew install bedtools parallel pigz r ruby

rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ $PWD/bin/kent-tools/

wget -O src/bowtie-1.1.2-linux-x86_64.zip \
     https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip/download \
    && unzip src/bowtie-1.1.2-linux-x86_64.zip -d $PWD/bin

wget -O src/tophat-2.1.1.Linux_x86_64.tar.gz \
     https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz \
    && tar zxvf src/tophat-2.1.1.Linux_x86_64.tar.gz -C $PWD/bin \
    && ln -s $PWD/bin/tophat-2.1.1.Linux_x86_64/samtools_0.1.18 $PWD/bin/tophat-2.1.1.Linux_x86_64/samtools

wget -O src/cufflinks-2.2.1.Linux_x86_64.tar.gz \
     http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz \
    && tar zxvf src/cufflinks-2.2.1.Linux_x86_64.tar.gz -C $PWD/bin

bundle

R --vanilla --quiet <<EOF
install.packages("BiocManager", repos="http://cran.r-project.org", Ncpus=`gnproc`)
BiocManager::install(c('remotes', 'yaml', 'impute', 'RColorBrewer', 'beeswarm'))
library(remotes)
install_github('renozao/pkgmaker@develop', force=T)
install_github('renozao/NMF@devel', force=T)
install_github('shka/samr', ref='test_multblock', force=T)
install_github('shka/R-SAMstrt', dependencies=F, upgrade='never', force=T)
EOF

gcc -o bin/_step1b_fastq2oneLine bin/_step1b_fastq2oneLine.c
gcc -o bin/_step1b_fastqs2oneLine bin/_step1b_fastqs2oneLine.c
gcc -o bin/_step1b_trimWithQCFilter bin/_step1b_trimWithQCFilter.c
