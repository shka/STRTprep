#!/usr/bin/env bash

bcl=$1

set -e

for lane in {1..4}; do
    java -jar ./bin/picard.jar CheckIlluminaDirectory \
         READ_STRUCTURE=85T6B \
         BASECALLS_DIR=./src/${bcl}/Data/Intensities/BaseCalls/ \
         LANES=${lane} > ./src/picard.CheckIlluminaDirectory.${bcl}.${lane} 2>&1
done

for lane in {1..4}; do
    java -Xmx7g -jar ./src/picard.jar IlluminaBasecallsToFastq \
         NUM_PROCESSORS=1 \
         READ_STRUCTURE=85T6B \
         BASECALLS_DIR=./src/${bcl}/Data/Intensities/BaseCalls/ \
         LANE=${lane} \
         OUTPUT_PREFIX=./src/${bcl}.lane${lane} \
         RUN_BARCODE=${bcl} \
         READ_NAME_FORMAT=ILLUMINA \
         INCLUDE_NON_PF_READS=false \
         COMPRESS_OUTPUTS=true \
         > ./src/picard.IlluminaBasecallsToFastq.${bcl}.${lane} 2>&1 &
done

wait
