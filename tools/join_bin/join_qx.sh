#!/bin/bash

# Helps to merge partitioned binary outputs in one file
# Author: Mario Morales-Hernandez

OUTPUT_DIR="${1:-../../outputParallel/bin}"
MERGE_DIR="${2:-../../outputParallel/merge}"
WHATMAT=QX
PATTERN="${3:-${WHATMAT}_[0-9][0-9]*_[0-9][0-9]*\.*}"
PATTERN0="${3:-${WHATMAT}_[0-9][0-9]*_00\.*}"

mkdir $MERGE_DIR

MERGE_DIR_BIN=$(echo "${MERGE_DIR}/bin")

make clean && make
./headerBinary "$OUTPUT_DIR/" "$MERGE_DIR_BIN"

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
	cat "${MERGE_DIR_BIN}/header.bin" >  "${MERGE_DIR_BIN}/${BASE_NAME}_00.out"
	while [ $p -lt $NP ]; do
		printf -v pp "%02d" "$p"
		tail -c +17 "${OUTPUT_DIR}/${BASE_NAME}_${pp}.out" > "${MERGE_DIR_BIN}/new.trash"
		cat "${MERGE_DIR_BIN}/new.trash" >> "${MERGE_DIR_BIN}/${BASE_NAME}_00.out"
		p=$((p+1))
	done
	p=0

done
rm "${MERGE_DIR_BIN}/new.trash"
rm "${MERGE_DIR_BIN}/header.bin"

exit 0
