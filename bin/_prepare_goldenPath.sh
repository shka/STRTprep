#!/usr/bin/env bash

ver=$1
tgz=src/${ver}_chromFa.tar.gz

curl -o $tgz "http://hgdownload.cse.ucsc.edu/goldenPath/$ver/bigZips/chromFa.tar.gz"

mkdir -p tmp
unpigz -c $tgz | tar -xOf - --exclude "*_*" | gawk 'BEGIN{p=""} /^>/{ print p $1 } !/^>/{ printf $1; p="\n" } END{ print }' | gfold -w 50 | pigz -c > tmp/$ver.fa.gz
