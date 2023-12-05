To build:
	module load gcc
	make clean && make

Use to convert a binary rmap file (rmap)
For example:
	./convert_rmap_file.sh

bin2ascii takes 3 arguments
	./bin2ascii_rmap $INPUT_DIR $OUTPUT_DIR $TYPE
	
$INPUT_DIR = Input file/folder directory of rmap binary files
$OUTPUT_DIR = Output file/folder directory of rmap ascii files
$TYPE = 1
