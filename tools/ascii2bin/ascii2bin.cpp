/**
 *
 * Helps to convert ascii file into binary format
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

std::pair<long, long> get_dims_2d(string& filepath)
{
	long num_cols = 0, num_rows = 0;
	bool first = false;
	ifstream infile(filepath.c_str());
		
	string line;

	while (getline(infile, line))
	{
		std::vector<std::string> row = split(line, ' ');

		if (!first)
		{
			num_cols = row.size();
			first = true;
		}
		num_rows++;
	}

	infile.close();

	return std::pair<long, long>(num_cols, num_rows);
}

void ascii2bin(string input_file, string output_file){
	std::tuple<long, long> dims = get_dims_2d(input_file);
	long nrows = std::get<1>(dims);
	long ncols = std::get<0>(dims);

    value_t *arr = new value_t [nrows*ncols];
	
	ifstream input(input_file.c_str());
	long i = 0;

	string line;

	while (input.good())
	{
		long j = 0;
		std::getline(input, line);
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
	input.close();

    ofstream output(output_file.c_str(), std::ios::binary);
	value_t put_rows_value = (value_t)(nrows);
	value_t put_cols_value = (value_t)(ncols);
			
	output.write((char*) &put_rows_value, sizeof(value_t));
	output.write((char*) &put_cols_value, sizeof(value_t));

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
		ascii2bin(input_dir, output_dir);
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
					
					ascii2bin(input_file, output_file);
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
