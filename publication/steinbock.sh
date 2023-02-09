#!/usr/bin/env bash

# change directory
BASEDIR=$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
cd ${BASEDIR}

# setup steinbock alias
shopt -s expand_aliases
alias steinbock="docker run -v ${BASEDIR}/data/steinbock:/data -u $(id -u):$(id -g) ghcr.io/bodenmillergroup/steinbock:0.16.0"

# image pre-processing
{ time steinbock preprocess imc images --hpf 50; } 2> steinbock_timing.txt

# image segmentation
{ time steinbock segment deepcell --minmax; } 2>> steinbock_timing.txt

# intensity measurement
{ time steinbock measure intensities; } 2>> steinbock_timing.txt

# regionprops measurement
{ time steinbock measure regionprops; } 2>> steinbock_timing.txt

# spatial cell graph construction
{ time steinbock measure neighbors --type expansion --dmax 4; } 2>> steinbock_timing.txt