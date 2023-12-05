#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/input/initqx/asc/case4.initqx"
OUTPUT_DIR="${PROJECT_ROOT}/input/initqx/bin/case4.initqx"
TYPE=1

./ascii2bin $INPUT_DIR $OUTPUT_DIR $TYPE

cd $START
