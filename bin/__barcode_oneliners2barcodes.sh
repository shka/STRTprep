#!/usr/bin/env bash

set -e
N=$(head -n 1 $1 | wc -w)
for i in `seq 1 $N`; do
    head -n 2 $1 | gawk '{if ($'$i' ~ /[0-9]/) print $'$i'; else print $'$i' "GGG"}' | tr '\n' '\t'
    echo
done
