#!/usr/bin/env bash
. bin/setup.sh
bundle
mkdir -p var/lib/R
R --vanilla --quiet <<EOF
install.packages(c('maptools', 'pvclust', 'Rmpi', 'snow'), repos='http://cran.us.r-project.org')
EOF
