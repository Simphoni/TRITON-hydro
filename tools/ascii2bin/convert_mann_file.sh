#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/input/mann/asc/case4.mann"
OUTPUT_DIR="${PROJECT_ROOT}/input/mann/bin/case4.mann"
TYPE=1

./ascii2bin $INPUT_DIR $OUTPUT_DIR $TYPE

cd $START