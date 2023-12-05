#!/bin/bash

# Helps to merge partitioned binary outputs in one file on the fly (while the simulation is running)
# Author: Mario Morales-Hernandez

OUTPUT_DIR="${1:-../../output/bin}"
MERGE_DIR="${2:-../../output/OTG}"
NTIMES="${3:-120}"

PATTERN="H_01"
BIN2ASCII_DIR="../bin2ascii/"
CURRENT_DIR=$PWD
echo $CURRENT_DIR

mkdir -p "$MERGE_DIR"

cd $BIN2ASCII_DIR
make clean && make
cp bin2ascii "$CURRENT_DIR/"
cd $CURRENT_DIR

while [ ! -f "${OUTPUT_DIR}/H_01_00.out" ]; do
	echo "Waiting for the first file..."
	sleep 10
done
#give time to write all the processes just in case
sleep 2 

# derive number of processes
NP=0
FILES=$(ls -p  $OUTPUT_DIR | grep -v / | grep -E "${PATTERN}")
for filename in $FILES; do
    x=$(echo $filename | cut -d'_' -f 3 | cut -d'.' -f 1)
    if [ $x -gt $NP ]; then
    	NP=$x
    fi
done
NP=$((10#$NP + 1))

make clean && make
./headerBinary "$OUTPUT_DIR/" "$MERGE_DIR/bin/"


t=1
while [ $t -le $NTIMES ]; do
	printf -v tt "%02d" "$t"
	BASE_NAME="H_$tt"
	while [ ! -f "${OUTPUT_DIR}/${BASE_NAME}_00.out" ]; do
		echo "Waiting for file" "${OUTPUT_DIR}/${BASE_NAME}_00.out ..."
		sleep 10
	done
	#give time to write all the processes and variables (just in case)
	sleep 2 

	for i in {1..3}; do
		case $i in
			1)
			BASE_NAME="H_${tt}"
			;;

			2)
			BASE_NAME="QX_${tt}"
			;;
			
			3)
			BASE_NAME="QY_${tt}"
			;;
		esac

		p=0
		cat "${MERGE_DIR}/bin/header.bin" >  "${MERGE_DIR}/bin/${BASE_NAME}_00.out"
		while [ $p -lt $NP ]; do
			printf -v pp "%02d" "$p"
			#remove the 16 first bytes that contain the information about the number of rows and columns in the binary file
			tail -c +17 "${OUTPUT_DIR}/${BASE_NAME}_${pp}.out" > "${MERGE_DIR}/bin/new.temp"
			cat "${MERGE_DIR}/bin/new.temp" >> "${MERGE_DIR}/bin/${BASE_NAME}_00.out"
			p=$((p+1))

		done
		./bin2ascii "${MERGE_DIR}/bin/${BASE_NAME}_00.out" "${MERGE_DIR}/asc/${BASE_NAME}_00.out" 1
		echo -e "[\e[32mOK\e[0m]" "File" "${MERGE_DIR}/bin/${BASE_NAME}_00.out" "converted"
	done
	t=$((t+1))

done
rm "${MERGE_DIR}/bin/header.bin"
rm "${MERGE_DIR}/bin/new.temp"
