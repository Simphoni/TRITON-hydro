To build:
	module load gcc
	make clean && make

Use to convert a binary dem file or all binary dem files in a directory
For example:
	./convert_dem_file.sh
	./convert_dem_folder.sh

bin2ascii_dem takes 3 arguments
	./bin2ascii_dem $INPUT_DIR $OUTPUT_DIR $TYPE
	
$INPUT_DIR = Input file/folder directory of Binary dem files
$OUTPUT_DIR = Output file/folder directory of Ascii dem files
$TYPE = 1 for file, 2 for folder