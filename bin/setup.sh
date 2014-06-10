#!/usr/bin/env bash

module load bioinfo-tools
module load samtools/0.1.19
module load bowtie/1.0.0
module load tophat/2.0.10
module load R/3.0.1
module load openmpi/1.7
module load ruby/2.1.0
module load BEDTools/2.17.0
module load ucsc/140527

export GEM_HOME=$PWD/var/lib/gems
export RUBYLIB=$GEM_HOME
export PATH=$GEM_HOME/bin:$PATH
export R_LIBS=$PWD/var/lib/R
