#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/output/bin/QX_01_00.out"
FILE_TYPE=2
TYPE=1

./discharge2velocity $INPUT_DIR $FILE_TYPE $TYPE

cd $START
