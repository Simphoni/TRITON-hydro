import os
import sys
import struct
import numpy as np
from osgeo import gdal

HEADER = 2
def bin2geotiff(input_file, output_file, xll, yll, cellsize):
    with open(input_file, 'rb') as input:
        dim_arr = struct.unpack('d' * HEADER, input.read(HEADER * struct.calcsize('d')))
        nrows = int(dim_arr[0])
        ncols = int(dim_arr[1])
        #arr = struct.unpack('d' * nrows * ncols, input.read(struct.calcsize('d') * nrows * ncols))
        arr = np.fromfile(input, dtype=np.float64, count=nrows*ncols).reshape((nrows, ncols))

    print(f'Input file read')
    print(type(arr), ncols, nrows)
    # Create a new GeoTIFF file with a single band
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(output_file, ncols, nrows, 1, gdal.GDT_Float64)

    # Set the georeferencing information (if available)
    geotransform = (xll, cellsize, 0.0, yll, 0.0, -cellsize)
    dataset.SetGeoTransform(geotransform)

    # Write the data to the raster band
    band = dataset.GetRasterBand(1)
    band.WriteArray(arr, 0, 0)

    # Clean up
    dataset = None

def main():

    # Check if 5 arguments are provided
    if len(sys.argv) != 6:
        print("Error: Please provide input_file, output_file, XLL, YLL, and CELLSIZE as arguments")
        sys.exit(1)

    # Get input and output file paths from command line arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Get XLL, YLL, and cellsize from command line arguments
    xll = float(sys.argv[3])
    yll = float(sys.argv[4])
    cellsize = float(sys.argv[5])

    # Create output directory if it doesn't exist
    parent_dir = os.path.dirname(output_file)
    if not os.path.exists(parent_dir):
        os.makedirs(parent_dir)

    print(f'Converting {input_file} to {output_file}...')

    # Convert binary file to GeoTIFF
    bin2geotiff(input_file, output_file,xll,yll,cellsize)

    print('Done.')


if __name__ == '__main__':
    main()

