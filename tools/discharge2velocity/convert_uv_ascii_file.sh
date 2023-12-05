#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/output/asc/QX_01_00.out"
FILE_TYPE=1
TYPE=1

./discharge2velocity $INPUT_DIR $FILE_TYPE $TYPE

cd $START
