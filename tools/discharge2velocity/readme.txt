To build:
	module load gcc
	make clean && make

Use to convert a discharge file or all discharge files in a directory to velocity
For example:
	./convert_uv_ascii.file.sh
	./convert_uv_binary_folder.sh

discharge2velocity takes 3 arguments
	./discharge2velocity $INPUT_DIR $FILE_TYPE $TYPE
	
$INPUT_DIR = Input file/folder directory of discharge
$FILE_TYPE = 1 for Ascii, 2 for Binary
$TYPE = 1 for file, 2 for folder