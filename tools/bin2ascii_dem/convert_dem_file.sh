#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/input/dem/bin/case4.dem"
OUTPUT_DIR="${PROJECT_ROOT}/input/dem/asc/case4.dem"
TYPE=1

./bin2ascii_dem $INPUT_DIR $OUTPUT_DIR $TYPE

cd $START