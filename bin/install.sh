#!/usr/bin/env bash
. bin/setup.sh
gem install bundler
bundle
mkdir -p var/lib/R
## R CMD javareconf -e
export JAVA_HOME=/sw/comp/java/x86_64/sun_jdk1.7.0_25
export JAVA=$JAVA_HOME/jre/bin/java
export JAVAC=$JAVA_HOME/bin/javac
export JAVAH=$JAVA_HOME/bin/javah
export JAR=$JAVA_HOME/bin/jar
export JAVA_LIBS="-L$JAVA_HOME/jre/lib/amd64/server -ljvm"
export JAVA_CPPFLAGS="-I$JAVA_HOME/include -I$JAVA_HOME/include/linux"
export JAVA_LD_LIBRARY_PATH=$JAVA_HOME/jre/lib/amd64/server
R --vanilla --quiet <<EOF
install.packages(c('maptools', 'pvclust', 'Rmpi', 'snow', 'xlsx'), repos='http://cran.us.r-project.org')
EOF

