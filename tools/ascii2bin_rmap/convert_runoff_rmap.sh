#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/input/runoff/case03_runoff.rmap"
OUTPUT_DIR="${PROJECT_ROOT}/input/runoff/case03_runoff.rmap.bin"
TYPE=1

./ascii2bin_rmap $INPUT_DIR $OUTPUT_DIR $TYPE

cd $START
