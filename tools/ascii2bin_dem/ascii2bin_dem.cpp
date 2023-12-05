/**
 *
 * Helps to convert ascii dem file into binary format
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
#include <vector>
#include <sstream>

using namespace std;

typedef double value_t;
#define FILE 1
#define FOLDER 2
#define DEM_HEADER_SIZE 6

std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &elems)
{
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim))
	{
		elems.push_back(item);
	}
	return elems;
}

std::vector<std::string> split(const std::string &s, char delim)
{
	std::vector<std::string> elems;
	return split(s, delim, elems);
}


void ascii2bin_dem(string input_file, string output_file){
	ifstream ifs(input_file.c_str());
	long line_num = 0;
	
	long nrows, ncols;
	value_t xllcorner, yllcorner, cellsize, no_data_value;

	while (true)
	{
		string line;
		getline(ifs, line);
		if (!ifs)
			break;

		line_num++;
		if (line_num > DEM_HEADER_SIZE)
			break;

		std::vector<std::string> tokens = split(line, ' ');
		const char* value = (*(tokens.end() - 1)).c_str();

		if (!line.empty() && line[0] != '#')
		{
			switch (line_num)
			{
				case 1:
				{
					ncols = atoi(value);
					break;
				}
				case 2:
				{
					nrows = atoi(value);
					break;
				}
				case 3:
				{
					xllcorner = atof(value);
					break;
				}
				case 4:
				{
					yllcorner = atof(value);
					break;
				}
				case 5:
				{
					cellsize = atof(value);
					break;
				}
				case 6:
				{
					no_data_value = atof(value);
					break;
				}
				default:
				{

				}
			}
		}
	}
	ifs.close();
	
	value_t *arr = new value_t [nrows*ncols];
	long i = 0;
	ifstream infile(input_file.c_str());
	std::string line;

	long line_number = 0;

	while (infile.good())
	{
		std::getline(infile, line);

		line_number++;
		if (line_number <= DEM_HEADER_SIZE)
			continue;

		long j = 0;
		std::vector<std::string> row = split(line, ' ');
		std::string val;
		std::vector<std::string>::iterator strit = row.begin();

		for (; strit != row.end(); strit++, j++)
		{
			val = *strit;

			arr[(ncols * i) + j] = (val.find(".") != std::string::npos) ? (value_t)atof(val.c_str()) : (value_t)atoi(val.c_str());
		}
		i++;
	}
	infile.close();

    ofstream output(output_file.c_str(), std::ios::binary);
	value_t put_rows_value = (value_t)(nrows);
	value_t put_cols_value = (value_t)(ncols);
	value_t put_no_data_value = (value_t)(no_data_value);

	output.write((char*) &put_cols_value, sizeof(value_t));
	output.write((char*) &put_rows_value, sizeof(value_t));
	output.write((char*) &xllcorner, sizeof(value_t));
	output.write((char*) &yllcorner, sizeof(value_t));
	output.write((char*) &cellsize, sizeof(value_t));
	output.write((char*) &put_no_data_value, sizeof(value_t));

	output.write((char*)&arr[0], nrows*ncols * sizeof(value_t));
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
		ascii2bin_dem(input_dir, output_dir);
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
					
					ascii2bin_dem(input_file, output_file);
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
