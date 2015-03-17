#!/usr/bin/env bash

. bin/setup.sh

brew list --versions

echo "----"

bundle list

echo "----"

R --vanilla --quiet <<EOF
sapply(library()[['results']][, 'Package'],
       function(p) { t <- packageVersion(p); paste(t, collapse='.') })
EOF
