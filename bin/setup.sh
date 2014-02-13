#!/usr/bin/env bash
module load bowtie/1.0.0
module load ruby/2.1.0-preview1
module load samtools/0.1.19
module load tophat/2.0.9

export GEM_HOME=$PWD/var/lib/gems
export RUBYLIB=$GEM_HOME
export PATH=$GEM_HOME/bin:$PATH
export TMPDIR=$HOME/opt/tmp
