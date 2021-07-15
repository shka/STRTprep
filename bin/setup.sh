#!/usr/bin/env bash

eval $($PWD/homebrew/bin/brew shellenv)

export LC_ALL=C
export PATH=$PWD/bin:$PWD/homebrew/opt/ruby@2.7/bin:$PWD/homebrew/lib/ruby/gems/2.7.0/bin:$PWD/bin/bowtie-1.1.2:$PWD/bin/tophat-2.1.1.Linux_x86_64:$PWD/bin/cufflinks-2.2.1.Linux_x86_64:$PWD/bin/kent-tools:$PATH
