#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/output/bin/"
OUTPUT_DIR="${PROJECT_ROOT}/output/asc/"
TYPE=2

./bin2ascii $INPUT_DIR $OUTPUT_DIR $TYPE

cd $START