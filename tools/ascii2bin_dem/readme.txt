To build:
	module load gcc
	make clean && make

Use to convert a ascii dem file or all ascii dem files in a directory
For example:
	./convert_dem_file.sh
	./convert_dem_folder.sh

ascii2bin_dem takes 3 arguments
	./ascii2bin_dem $INPUT_DIR $OUTPUT_DIR $TYPE
	
$INPUT_DIR = Input file/folder directory of Ascii dem files
$OUTPUT_DIR = Output file/folder directory of Binary dem files
$TYPE = 1 for file, 2 for folder