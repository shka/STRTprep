#!/usr/bin/env bash

export LC_ALL=C
export GEM_HOME=$PWD/vendor/gems
export RUBYLIB=$GEM_HOME
export R_LIBS=$PWD/vendor/Rlibs
export PATH=$GEM_HOME/bin:.homebrew/bin:$PATH

mkdir -p $R_LIBS $GEM_HOME

module use ~katay/.modulefiles
module load bioinfo-tools
module load ruby/2.1.0
module load gnuparallel/20140222
module load BEDTools/2.21.0
module load samtools/0.1.19
module load bowtie/1.1.0
module load tophat/2.0.12
module load ucsc-utilities/v287
module load R/3.2.3
module load cufflinks/2.1.1
