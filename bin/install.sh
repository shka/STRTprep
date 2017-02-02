#!/usr/bin/env bash

. bin/setup.sh

mkdir -p $R_LIBS $GEM_HOME

gem install bundler
bundle
