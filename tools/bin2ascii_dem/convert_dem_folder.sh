#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/input/dem/bin/"
OUTPUT_DIR="${PROJECT_ROOT}/input/dem/asc/"
TYPE=1

./bin2ascii_dem $INPUT_DIR $OUTPUT_DIR $TYPE

cd $START