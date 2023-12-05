#!/bin/bash

# Helps to merge partitioned ascii outputs in one file
# Author: Ryan Marshall

OUTPUT_DIR="${1:-../../output/asc}"
WHATMAT=H
PATTERN="${2:-${WHATMAT}_[0-9][0-9]*_[0-9][0-9]*\.*}"
PATTERN0="${2:-${WHATMAT}_[0-9][0-9]*_00\.*}"
DELETE_OLD=${3:-1}

FILES=$(ls -p  $OUTPUT_DIR | grep -v / | grep -E "${PATTERN}")
PROC0_FILES=$(ls -p  $OUTPUT_DIR | grep -v / | grep -E "${PATTERN0}")
NP=0

# derive number of processes
for filename in $FILES; do
    x=$(echo $filename | cut -d'_' -f 3 | cut -d'.' -f 1)
    if [ $x -gt $NP ]; then
    	NP=$x
    fi
done

# convert MPI rank id (0-based) to number of processes (1-based)

NP=$((10#$NP + 1))

p=0

# only loop through the filenames for process 0.  We already know NP, 
#  so we can guess the remaining filenames
for filename in $PROC0_FILES; do
	BASE_NAME=$(echo $filename | cut -d'_' -f 1,2)
	echo "$BASE_NAME"
	while [ $p -lt $NP ]; do
		printf -v pp "%02d" "$p"
		cat "${OUTPUT_DIR}/${BASE_NAME}_${pp}.out" >> "${OUTPUT_DIR}/${BASE_NAME}.out"
		if [ $DELETE_OLD -gt 0 ]; then
			rm "${OUTPUT_DIR}/${BASE_NAME}_${pp}.out"
		fi
		p=$((p+1))
	done
	p=0

done

exit 0
