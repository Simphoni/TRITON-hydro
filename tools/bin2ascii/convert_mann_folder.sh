#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/input/mann/bin/"
OUTPUT_DIR="${PROJECT_ROOT}/input/mann/asc/"
TYPE=2

./bin2ascii $INPUT_DIR $OUTPUT_DIR $TYPE

cd $START