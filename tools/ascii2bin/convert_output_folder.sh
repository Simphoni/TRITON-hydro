#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/output/asc/"
OUTPUT_DIR="${PROJECT_ROOT}/output/bin/"
TYPE=2

./ascii2bin $INPUT_DIR $OUTPUT_DIR $TYPE

cd $START