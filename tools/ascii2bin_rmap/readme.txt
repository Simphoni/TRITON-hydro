To build:
	module load gcc
	make clean && make

Use to convert an rmap file (*.rmap)
For example:
	./convert_runoff_rmap.sh

ascii2bin_rmap takes 3 arguments
	./ascii2bin $INPUT_DIR $OUTPUT_DIR $TYPE
	
$INPUT_DIR = Input file/folder directory of rmap ascii file
$OUTPUT_DIR = Output file/folder directory of rmap binary file
$TYPE = 1
