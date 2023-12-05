#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/input/initqy/bin/case4.initqy"
OUTPUT_DIR="${PROJECT_ROOT}/input/initqy/asc/case4.initqy"
TYPE=1

./bin2ascii $INPUT_DIR $OUTPUT_DIR $TYPE

cd $START
