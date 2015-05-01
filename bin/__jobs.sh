#!/usr/bin/env bash

set -e

. bin/setup.sh
rake -m -j `gnproc` qc > qc.log 2>&1
rake gene > gene.log 2>&1
rake > tfe.log 2>&1

