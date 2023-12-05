#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <dirent.h>
#include <cstring>
#include <sys/stat.h>

using namespace std;

typedef double value_t;
#define HEADER 2


int main(int argc, char* argv[])
{			
	string output_dir(argv[1]);
	int ncols, nrows;

	nrows=atoi(argv[2]);
	ncols=atoi(argv[3]);

	DIR* dir_ = opendir(output_dir.empty() ? "." : output_dir.c_str());
	if( !dir_) {
		mkdir(output_dir.c_str(), S_IRWXU);
	}
	else closedir(dir_);

	string output_file = output_dir + "/header.bin";

	ofstream output(output_file.c_str(), std::ios::binary);
	value_t put_rows_value = (value_t)(nrows);
	value_t put_cols_value = (value_t)(ncols);
	
	output.write((char*) &put_rows_value, sizeof(value_t));
	output.write((char*) &put_cols_value, sizeof(value_t));

   output.close();


    return 0;
}
