To build:
	module load gcc
	make clean && make

Use to convert a ascii file or all ascii files in a directory
For example:
	./convert_output_file.sh
	./convert_output_folder.sh

ascii2bin takes 3 arguments
	./ascii2bin $INPUT_DIR $OUTPUT_DIR $TYPE
	
$INPUT_DIR = Input file/folder directory of Ascii files
$OUTPUT_DIR = Output file/folder directory of Binary files
$TYPE = 1 for file, 2 for folder