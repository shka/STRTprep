#!/usr/bin/env bash

set -e

if ! [ -d '.homebrew' ]; then
    if [ `uname` = 'Linux' ]; then
	git clone https://github.com/Homebrew/linuxbrew.git .homebrew
    else
	mkdir .homebrew
	curl -L https://github.com/Homebrew/homebrew/tarball/master | tar xz --strip 1 -C .homebrew
    fi
fi

. bin/setup.sh

brew update

brew tap homebrew/science

# bedtools, not later than 2.22.0
rm .homebrew/Library/Taps/homebrew/homebrew-science/bedtools.rb
git --work-tree .homebrew/Library/Taps/homebrew/homebrew-science --git-dir .homebrew/Library/Taps/homebrew/homebrew-science/.git checkout 5392a9e bedtools.rb
brew install bedtools

brew install ruby
brew install coreutils
brew install fastq-tools
brew install parallel
brew install pigz
brew install gawk
brew install samtools
brew install bowtie
brew install tophat --without-bowtie2
# brew install edirect
brew install https://raw.githubusercontent.com/Homebrew/homebrew-science/fbf8b1f20c27baa29c24431a03cf30868d6cc933/kent-tools.rb
brew install R --with-openblas --without-tcltk --without-x11

R CMD javareconf

gem install bundler
bundle

R --vanilla --quiet <<EOF
source("http://bioconductor.org/biocLite.R")
biocLite(c('Heatplus', 'xlsx', 'devtools', 'samr'), ask=F)
library(devtools)
install_github('shka/samr', ref='test_multblock')
install_github('shka/R-SAMstrt')
EOF
