#!/usr/bin/env bash

set -e

ver=$1
tgz=src/${ver}_chromFa.tar.gz

if [ $ver = 'hg38as' ]; then
    curl -o $tgz "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.chroms.tar.gz"
elif [ $ver = 'hg38' ]; then
    curl -o $tgz "http://hgdownload.cse.ucsc.edu/goldenPath/$ver/bigZips/hg38.chromFa.tar.gz"
else
    curl -o $tgz "http://hgdownload.cse.ucsc.edu/goldenPath/$ver/bigZips/chromFa.tar.gz"
fi

mkdir -p tmp
unpigz -c $tgz | tar -xOf - --exclude "*_*" | gawk 'BEGIN{p=""} /^>/{ print p $1 } !/^>/{ printf $1; p="\n" } END{ print }' | gfold -w 50 | pigz -c > tmp/$ver.fa.gz
