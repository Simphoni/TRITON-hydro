#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/output/bin/H_01_00.out"
OUTPUT_DIR="${PROJECT_ROOT}/output/asc/H_01_00.out"
TYPE=1

./bin2ascii $INPUT_DIR $OUTPUT_DIR $TYPE

cd $START