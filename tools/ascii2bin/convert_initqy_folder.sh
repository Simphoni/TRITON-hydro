#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/input/initqy/asc/"
OUTPUT_DIR="${PROJECT_ROOT}/input/initqy/bin/"
TYPE=2

./ascii2bin $INPUT_DIR $OUTPUT_DIR $TYPE

cd $START
