#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/input/inith/bin/case4.inith"
OUTPUT_DIR="${PROJECT_ROOT}/input/inith/asc/case4.inith"
TYPE=1

./bin2ascii $INPUT_DIR $OUTPUT_DIR $TYPE

cd $START