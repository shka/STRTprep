#!/usr/bin/env bash

. bin/setup.sh

mkdir -p .homebrew
curl -L https://github.com/Homebrew/homebrew/tarball/master | tar xz --strip 1 -C .homebrew

brew tap homebrew/science

brew install ruby samtools bowtie bedtools openmpi edirect
brew install https://raw.githubusercontent.com/Homebrew/homebrew-science/fbf8b1f20c27baa29c24431a03cf30868d6cc933/kent-tools.rb
brew install tophat --without-bowtie2
brew install R --with-openblas --without-tcltk --without-x11

R CMD javareconf

gem install bundler
bundle

R --vanilla --quiet <<EOF
install.packages('Rmpi', repos='http://cran.r-project.org', configure.args=sprintf("--with-mpi=%s/.homebrew --with-Rmpi-type=MPICH2", getwd()))
install.packages(c('maptools', 'pvclust', 'xlsx', 'snow'), repos='http://cran.r-project.org')
EOF

