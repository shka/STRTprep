#!/usr/bin/env bash

module load bioinfo-tools
module load samtools/0.1.19
module load bowtie/1.0.0
module load tophat/2.0.10
module load ruby/2.1.0
module load bedtools/2.17.0
module load ucsc

export GEM_HOME=$PWD/var/lib/gems
export RUBYLIB=$GEM_HOME
export PATH=$GEM_HOME/bin:$PATH
