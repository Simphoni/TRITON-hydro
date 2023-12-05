#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/output/asc/"
FILE_TYPE=1
TYPE=2

./discharge2velocity $INPUT_DIR $FILE_TYPE $TYPE

cd $START