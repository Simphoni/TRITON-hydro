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
#include <gdal/gdal.h>
#include <gdal/gdal_priv.h>

#define xll 0.0
#define yll 0.0
#define cellsize 30.0

using namespace std;

typedef double value_t;
#define HEADER 2

void bin2geotiff(string input_file, string output_file){
	ifstream input(input_file.c_str(), ios::binary);

	value_t *dim_arr = new value_t[HEADER];
	input.read((char*) dim_arr, HEADER * sizeof(value_t)); 

	long nrows = (long) dim_arr[0];
	long ncols = (long) dim_arr[1];

    value_t *arr = new value_t [nrows*ncols];

	input.seekg(HEADER * sizeof(value_t), std::ios::beg);
    input.read((char*) arr, sizeof(value_t) * nrows * ncols);
    input.close();

// Register all the GDAL drivers
    GDALAllRegister();

    // Create a new GeoTIFF file with a single band
    GDALDriver *driver = GetGDALDriverManager()->GetDriverByName("GTiff");
    GDALDataset *dataset = driver->Create(output_file.c_str(), ncols, nrows, 1, GDT_Float64, NULL);

    // Set the georeferencing information (if available)
    double geotransform[6] = {xll, cellsize, 0.0, yll, 0.0, cellsize};
    dataset->SetGeoTransform(geotransform);

    // Write the data to the raster band
    GDALRasterBand *band = dataset->GetRasterBand(1);
    band->RasterIO(GF_Write, 0, 0, ncols, nrows, arr, ncols, nrows, GDT_Float64, 0, 0);

    // Clean up
    delete[] arr;
    delete[] dim_arr;
    GDALClose(dataset);
}

int main(int argc, char* argv[])
{			
	string input_dir(argv[1]);
	string output_dir(argv[2]);
	
	string parent_dir = output_dir.substr(0, output_dir.find_last_of("/"));
	DIR* dir_ = opendir(parent_dir.empty() ? "." : parent_dir.c_str());
	if( !dir_) {
		mkdir(parent_dir.c_str(), S_IRWXU);
	}
	else closedir(dir_);
	cout << "Converting " << input_dir << endl;
	bin2geotiff(input_dir, output_dir);
	
   return 0;
}
