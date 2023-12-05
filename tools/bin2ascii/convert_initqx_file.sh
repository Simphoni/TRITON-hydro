#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/input/initqxu/bin/case4.initqxu"
OUTPUT_DIR="${PROJECT_ROOT}/input/initqxu/asc/case4.initqxu"
TYPE=1

./bin2ascii $INPUT_DIR $OUTPUT_DIR $TYPE

cd $START
