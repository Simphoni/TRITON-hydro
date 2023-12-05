/**
 *
 * Helps to convert discharge to velocity
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
#define ASCII 1
#define BINARY 2
#define HEXTRA 0.001
#define HEADER 2

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

std::pair<int, int> get_dims_2d(string& filepath)
{
	int num_cols = 0, num_rows = 0;
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

	return std::pair<int, int>(num_cols, num_rows);
}

void discharge2velocity_ascii(string input_file, string input_file_h, string output_file){
	std::tuple<int, int> dims = get_dims_2d(input_file);
	int rows = std::get<1>(dims);
	int cols = std::get<0>(dims);
	
	ifstream input(input_file.c_str());
	ifstream input_h(input_file_h.c_str());
	int i = 0;

	string line, line_h;

    ofstream output(output_file.c_str());

	while (input.good() && input_h.good())
	{
		int j = 0;
		std::getline(input, line);
		std::getline(input_h, line_h);
		std::vector<std::string> row = split(line, ' ');
		std::vector<std::string> row_h = split(line_h, ' ');
		std::string val;
		std::string val_h;
		std::vector<std::string>::iterator strit = row.begin();
		std::vector<std::string>::iterator strit_h = row_h.begin();

		for (; strit != row.end() && strit_h != row_h.end(); strit++, strit_h++, j++)
		{
			val = *strit;
			val_h = *strit_h;
			value_t value = (val.find(".") != std::string::npos) ? (value_t)atof(val.c_str()) : (value_t)atoi(val.c_str());
			value_t value_h = (val_h.find(".") != std::string::npos) ? (value_t)atof(val_h.c_str()) : (value_t)atoi(val_h.c_str());
			output << std::setprecision(6) << std::fixed;
			if(value_h < HEXTRA)
			{
				output << 0.0;
			}
			else
			{
				output << (value/value_h);
			}
			if(j < cols-1){
                output << " ";
            }
		}
		i++;
		output << endl;
	}
	input.close();
	input_h.close();
    output.close();
}

void discharge2velocity_binary(string input_file, string input_file_h, string output_file){
	ifstream input(input_file.c_str(), ios::binary);

	value_t *dim_arr = new value_t[HEADER];
	input.read((char*) dim_arr, HEADER * sizeof(value_t)); 

	int row = (int) dim_arr[0];
	int col = (int) dim_arr[1];

    value_t *arr = new value_t [row*col];

	input.seekg(HEADER * sizeof(value_t), std::ios::beg);
    input.read((char*) arr, sizeof(value_t) * row * col);
    input.close();
	
	ifstream input_h(input_file_h.c_str(), ios::binary);
	
	value_t *arr_h = new value_t [row*col];
	input_h.seekg(HEADER * sizeof(value_t), std::ios::beg);
    input_h.read((char*) arr_h, sizeof(value_t) * row * col);
    input_h.close();
	
	for(int i=0; i< row*col; i++)
	{
		if(arr_h[i] < HEXTRA)
		{
			arr[i] = 0.0;
		}
		else
		{
			arr[i] = arr[i]/arr_h[i];
		}
	}

    ofstream output(output_file.c_str(), std::ios::binary);
	value_t put_rows_value = (value_t)(row);
	value_t put_cols_value = (value_t)(col);
			
	output.write((char*) &put_rows_value, sizeof(value_t));
	output.write((char*) &put_cols_value, sizeof(value_t));

	output.write((char*)&arr[0], row*col * sizeof(value_t));
    output.close();

    delete[] arr;
	delete[] arr_h;
}

int main(int argc, char* argv[])
{			
	string input_dir(argv[1]);
	int file_type = atoi(argv[2]);
	int type = atoi(argv[3]);
	
	if (type == FILE)
	{
		cout << "Converting " << input_dir << endl;
		size_t index = 0;
		if(input_dir.find("QX_") != std::string::npos)
		{
			index = input_dir.find("QX_", index);
			string input_dir_h = input_dir.substr(0, index) + "H_" + input_dir.substr(index+3);
			string output_dir = input_dir.substr(0, index) + "U_" + input_dir.substr(index+3);
			
			if(file_type == ASCII){
				discharge2velocity_ascii(input_dir, input_dir_h, output_dir);
			}
			else
			{
				discharge2velocity_binary(input_dir, input_dir_h, output_dir);
			}
		}
		else if(input_dir.find("QY_") != std::string::npos)
		{
			index = input_dir.find("QY_", index);
			string input_dir_h = input_dir.substr(0, index+1) + "H_" + input_dir.substr(index+3);
			string output_dir = input_dir.substr(0, index+1) + "V_" + input_dir.substr(index+3);
			if(file_type == ASCII){
				discharge2velocity_ascii(input_dir, input_dir_h, output_dir);
			}
			else
			{
				discharge2velocity_binary(input_dir, input_dir_h, output_dir);
			}
		}
	}
	else if (type == FOLDER)
	{
		DIR *dir;
		struct dirent *ent;
		if ((dir = opendir (input_dir.c_str())) != NULL) {
			while ((ent = readdir (dir)) != NULL) {
				if( strcmp(ent->d_name, ".") != 0 && strcmp(ent->d_name, "..") != 0 ){
					string str (ent->d_name);
					if (str.find("QX_") != std::string::npos || str.find("QY_") != std::string::npos) {
						cout << "Converting " << str << endl;
					
						string input_file = input_dir + str;
						string input_file_h = input_dir + "H_" + str.substr(3);
						string output_file = input_dir + "U_" + str.substr(3);
						if(str.find("QY_") != std::string::npos)
						{
							output_file = input_dir + "V_" + str.substr(3);
						}
					
						if(file_type == ASCII){
							discharge2velocity_ascii(input_file, input_file_h, output_file);
						}
						else
						{
							discharge2velocity_binary(input_file, input_file_h, output_file);
						}
					}
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
