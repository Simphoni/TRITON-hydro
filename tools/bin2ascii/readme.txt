To build:
	module load gcc
	make clean && make

Use to convert a binary file or all binary files in a directory
For example:
	./convert_output_file.sh
	./convert_output_folder.sh

bin2ascii takes 3 arguments
	./bin2ascii $INPUT_DIR $OUTPUT_DIR $TYPE
	
$INPUT_DIR = Input file/folder directory of Binary files
$OUTPUT_DIR = Output file/folder directory of Ascii files
$TYPE = 1 for file, 2 for folder