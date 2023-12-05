#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/input/runoff/case03_runoff.rmap.bin"
OUTPUT_DIR="${PROJECT_ROOT}/input/runoff/case03_runoff.rmap"
TYPE=1

./bin2ascii_rmap $INPUT_DIR $OUTPUT_DIR $TYPE

cd $START
