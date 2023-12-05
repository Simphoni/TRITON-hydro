/**
 *
 * Helps to convert binary file into ascii format
 * Author: Md Bulbul Sharif
 *
 **/

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

void readHeaderBinary(string input_file, int *row, int *col){
	ifstream input(input_file.c_str(), ios::binary);

	value_t *dim_arr = new value_t[HEADER];
	input.read((char*) dim_arr, HEADER * sizeof(value_t)); 

	*row = (int) dim_arr[0];
	*col = (int) dim_arr[1];

   input.close();
}

int main(int argc, char* argv[])
{			
	string input_dir(argv[1]);
	string output_dir(argv[2]);
	
	DIR* dir_ = opendir(output_dir.empty() ? "." : output_dir.c_str());
	if( !dir_) {
		mkdir(output_dir.c_str(), S_IRWXU);
	}
	else closedir(dir_);

	DIR *dir;
	struct dirent *ent;
	int totalrows=0,totalcols=0;
	int row,col;
	int count=0;
	if ((dir = opendir (input_dir.c_str())) != NULL) {
		while ((ent = readdir (dir)) != NULL) {
			if( strcmp(ent->d_name, ".") != 0 && strcmp(ent->d_name, "..") != 0 ){
				string str (ent->d_name);
				if(str.find("_01_") != string::npos){
					//cout << "Reading " << str << endl;
					string input_file = input_dir + str;
					readHeaderBinary(input_file, &row,&col);
					totalrows+=row;
					totalcols+=col;
					count++;
				}
			}
		}
		closedir (dir);
	} else {
		perror ("");
		return 1;
	}
	
	//printf("%d %d\n",totalrows/3,totalcols/count);

	string output_file = output_dir + "/header.bin";
	
	ofstream output(output_file.c_str(), std::ios::binary);
	//divided by three: H, QX, QY. 
	value_t put_rows_value = (value_t)(totalrows/3);
	//divided by count: row-wise decomposition
	value_t put_cols_value = (value_t)(totalcols/count);
			
	output.write((char*) &put_rows_value, sizeof(value_t));
	output.write((char*) &put_cols_value, sizeof(value_t));

   output.close();


    return 0;
}
