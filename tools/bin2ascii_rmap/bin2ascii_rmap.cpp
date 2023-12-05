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

typedef int value_t;
#define FILE 1
#define FOLDER 2
#define HEADER 2

void bin2ascii(string input_file, string output_file){
	ifstream input(input_file.c_str(), ios::binary);

	value_t *dim_arr = new value_t[HEADER];
	input.read((char*) dim_arr, HEADER * sizeof(value_t)); 

	long nrows = (long) dim_arr[0];
	long ncols = (long) dim_arr[1];

    value_t *arr = new value_t [nrows*ncols];

	input.seekg(HEADER * sizeof(value_t), std::ios::beg);
    input.read((char*) arr, sizeof(value_t) * nrows * ncols);
    input.close();

    ofstream output(output_file.c_str());
    output << std::setprecision(6) << std::fixed;

    for(long i=0; i<nrows; i++){
        for(long j=0; j<ncols; j++){
            output << arr[i*ncols+j];
            if(j<ncols-1){
                output << " ";
            }
        }
        output << endl;
    }
    output.close();

    delete[] arr;
}

int main(int argc, char* argv[])
{			
	string input_dir(argv[1]);
	string output_dir(argv[2]);
	int type = atoi(argv[3]);
	
	if (type == FILE)
	{
		string parent_dir = output_dir.substr(0, output_dir.find_last_of("/"));
		DIR* dir_ = opendir(parent_dir.empty() ? "." : parent_dir.c_str());
		if( !dir_) {
			mkdir(parent_dir.c_str(), S_IRWXU);
		}
		else closedir(dir_);
		cout << "Converting " << input_dir << endl;
		bin2ascii(input_dir, output_dir);
	}
	else if (type == FOLDER)
	{
		DIR* dir_ = opendir(output_dir.empty() ? "." : output_dir.c_str());
		if( !dir_) {
			mkdir(output_dir.c_str(), S_IRWXU);
		}
		else closedir(dir_);

		DIR *dir;
		struct dirent *ent;
		if ((dir = opendir (input_dir.c_str())) != NULL) {
			while ((ent = readdir (dir)) != NULL) {
				if( strcmp(ent->d_name, ".") != 0 && strcmp(ent->d_name, "..") != 0 ){
					string str (ent->d_name);
					cout << "Converting " << str << endl;
					
					string input_file = input_dir + str;
					string output_file = output_dir + str;
					
					bin2ascii(input_file, output_file);
				}
			}
			closedir (dir);
		} else {
			perror ("");
			return 1;
		}
	}

    return 0;
}
