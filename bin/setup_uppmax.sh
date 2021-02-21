#!/usr/bin/env bash
# [[file:~/Documents/org/Rackham.org::*STRTprep3%20configulation][STRTprep3 configulation:1]]
export LC_ALL=C
export GEM_HOME=$PWD/vendor/gems
export RUBYLIB=$GEM_HOME
export R_LIBS=$PWD/vendor/Rlibs
export PATH=$GEM_HOME/bin:/proj/uppstore2017139/private/.linuxbrew/bin:$PATH
export MANPATH="$(brew --prefix)/share/man:$MANPATH"
export INFOPATH="$(brew --prefix)/share/info:$INFOPATH"

mkdir -p $R_LIBS $GEM_HOME

module load git
module load ruby/2.5.0
module load gnuparallel
module load R/3.5.2
module load bioinfo-tools
module load BEDTools/2.21.0
module load samtools/0.1.19
module load bowtie/1.1.2
module load tophat/2.1.1
module load cufflinks/2.2.1
module load ucsc-utilities
# STRTprep3 configulation:1 ends here
