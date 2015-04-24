#!/usr/bin/env bash

set -e

mkdir -p bin/cufflinks-2.1.1
if ! [ -d '.homebrew' ]; then
    if [ `uname` = 'Linux' ]; then
	git clone https://github.com/Homebrew/linuxbrew.git .homebrew
	# cufflinks 2.1.1; 2.2.1 was unstable
	curl -s http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.1.1.Linux_x86_64.tar.gz | tar -C ./bin/cufflinks-2.1.1 --strip-component 1 -zxvf -
    else
	mkdir .homebrew
	curl -L https://github.com/Homebrew/homebrew/tarball/master | tar xz --strip 1 -C .homebrew
	# cufflinks 2.1.1; 2.2.1 was unstable
	curl -s http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.1.1.OSX_x86_64.tar.gz | tar -C ./bin/cufflinks-2.1.1 --strip-component 1 -zxvf -
    fi
fi

. bin/setup.sh

brew update
brew tap homebrew/science
brew install R --with-openblas --without-tcltk --without-x11

# bedtools, not later than 2.22.0
rm .homebrew/Library/Taps/homebrew/homebrew-science/bedtools.rb
git --work-tree .homebrew/Library/Taps/homebrew/homebrew-science --git-dir .homebrew/Library/Taps/homebrew/homebrew-science/.git checkout 5392a9e bedtools.rb
brew install bedtools

brew install ruby
brew install coreutils
brew install parallel
brew install pigz
brew install gawk
brew install samtools-0.1
brew install bowtie
brew install tophat --with-bowtie
brew install kent-tools
brew install md5sha1sum

gem install bundler
bundle

R --vanilla --quiet <<EOF
source("http://bioconductor.org/biocLite.R")
biocLite(c('Heatplus', 'devtools', 'samr'), ask=F)
library(devtools)
install_github('shka/samr', ref='test_multblock')
install_github('shka/R-SAMstrt')
EOF


