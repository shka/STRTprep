#!/usr/bin/env bash

export LC_ALL=C
export GEM_HOME=$PWD/vendor/gems
export RUBYLIB=$GEM_HOME
export R_LIBS=$PWD/vendor/Rlibs
export PATH=$GEM_HOME/bin:$PATH

if [ $(uname -n | grep -e '.uppmax.uu.se$') ]; then
    . bin/setup_uppmax.sh
fi

