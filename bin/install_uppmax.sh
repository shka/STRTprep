#!/usr/bin/env bash

set -e

if ! [ -d '.homebrew' ]; then
    git clone https://github.com/Homebrew/linuxbrew.git .homebrew
fi

. bin/setup_uppmax.sh

brew tap homebrew/science
brew update

brew install coreutils
brew install pigz

gem install bundler
bundle

R --vanilla --quiet <<EOF
source("http://bioconductor.org/biocLite.R")
biocLite(c('Heatplus', 'devtools', 'samr'), ask=F)
library(devtools)
install_github('shka/samr', ref='test_multblock')
install_github('shka/R-SAMstrt')
EOF
