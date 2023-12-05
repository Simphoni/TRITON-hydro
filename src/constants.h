/** @file constants.h
 *  @brief Header containing the Constants class
 *
 *  This contains the subroutines and eventually any 
 *  macros, constants, etc. needed for Constants class
 *
 *  @author Mario Morales Hernandez
 *  @author Md Bulbul Sharif
 *  @author Tigstu T. Dullo
 *  @author Sudershan Gangrade
 *  @author Alfred Kalyanapu
 *  @author Sheikh Ghafoor
 *  @author Shih-Chieh Kao
 *  @author Katherine J. Evans
 *  @bug No known bugs.
 */



#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Constants
{
	typedef std::pair<int, int> dims_t;    /**< Custom type to define dimension. The first number represents the rows and the second number is the columns. */
	typedef std::vector<std::string> string_vector;    /**< Custom vector type that represents a vector of strings. */
	typedef std::string::value_type char_t;    /**< Custom type used for string utility. */
	typedef std::vector<std::pair<int, int>> sources_list_t;	/**< Custom vector type that contains each cell's index pair. The first number is the row index and the second number is the column index. */
	typedef unsigned long long ull;    /**< Custom data type to hold large number. */
}

typedef double value_t;    /**< Data type to represent floating-point number. It can be double or float. */

#define MPI_DATA_TYPE MPI_DOUBLE    /**< Represents MPI floating-point number. It can be MPI_DOUBLE or MPI_FLOAT. */
#define MAX_VALUE DBL_MAX    /**< Maximum value of a floating-point number. It can be DBL_MAX or FLT_MAX. */

#define INPUT_DIR "input"    /**< Deafult folder name containing all input files. */
#define OUTPUT_DIR "output"    /**< Deafult folder name containing all output files. */
#define CFG_DIR "cfg"    /**< Default folder name containing all configuration (cfg) files. */
#define BIN_DIR "bin"    /**< Default folder name containing binary files. */
#define ASCII_DIR "asc"    /**< Deafult folder name containing ascii files. */
#define TIME_SERIES_DIR "series"    /**< Deafult folder name containing time series outputs. */
#define DEFAULT_CFG "case4.cfg"    /**< Deafult configuration (cfg) file name. */

#define GHOST_CELL_PADDING 1    /**< Number of extra row and column to use besides each domain border. */
#define USE_MATRIX 0    /**< Use the whole matrix when performing MPI halo exchange. */
#define USE_HALO 1    /**< Use only halo rows bundle when performing MPI halo exchange. */
#define SRC_LOCATION 0    /**< Define to use flow locations. */
#define OBSERVATION_LOCATION 1    /**< Define to use observation cells. */

#define DEM_NCOLS_LINE 1    /**< Line 1 in DEM file represents number of columns. */
#define DEM_NROWS_LINE 2    /**< Line 2 in DEM file represents number of rows. */
#define DEM_XLL_CORNER_LINE 3    /**< Line 3 in DEM file represents X coordinate of the origin (by center or lower left corner of the cell). */
#define DEM_YLL_CORNER_LINE 4    /**< Line 4 in DEM file represents Y coordinate of the origin (by center or lower left corner of the cell). */
#define DEM_CELL_SIZE_LINE 5    /**< Line 5 in DEM file represents cell size. */
#define DEM_NODATA_VALUE_LINE 6    /**< Line 6 in DEM file represents the input values to be NoData in the output raster. */
#define DEM_HEADER_SIZE 6    /**< Number of headers in a DEM input file. */

#define BIN_ROW_ID 0    /**< First number or index 0 in a binary output file represents number of rows. */
#define BIN_COL_ID 1    /**< Second number or index 1 in a binary output file represents number of columns. */
#define BIN_DEFAULT_HEADER_SIZE 2    /**< Number of headers in a binary output file. */

#define H 0    /**< Water depth array position in vector. */
#define QX 1    /**< Flux X array position in vector. */
#define QY 2    /**< Flux Y array position in vector. */
#define NMAN 3    /**< Manning array position in vector. */
#define DEM 4    /**< DEM array position in vector. */
#define MAXH 5    /**< Max values of water depth array position in vector. */
#define RHSH0 6    /**< Partial water depth 1 array position in vector. */
#define RHSH1 7    /**< Partial water depth 2 array position in vector. */
#define RHSQX0 8    /**< Partial flux X 1 array position in vector. */
#define RHSQX1 9    /**< Partial flux X 2 array position in vector. */
#define RHSQY0 10    /**< Partial flux Y 1 array position in vector. */
#define RHSQY1 11    /**< Partial flux Y 2 array position in vector. */
#define SQRTH 12   /**< Square root of water depth array position in vector. */
#define HALO 13    /**< Halo cells array position in vector. */
#define DT 14    /**< Reduction values container array when calculating min time step size, position in vector. */
#define HYGT 15    /**< Time of flow values array position in vector. */
#define HYGV 16    /**< Flow values array position in vector. */
#define RUNIN 17    /**< Runoff intensity array position in vector. */
#define EXTBCV1 18    /**< External boundary condition's first variable array position in vector. */
#define EXTBCV2 19    /**< External boundary condition's second variable array position in vector. */

#define SRCP 0    /**< Flow locations index array position in vector. */
#define RUNID 1    /**< Runoff id array position in vector. */
#define BCRELATIVEINDEX 2    /**< Boundary cells index array after domain decomposition position in vector. */
#define BCTYPE 3    /**< Boundary condition cells type array position in vector. */
#define BCINDEXSTART 4    /**< Boundary condition's start index array position in vector. */
#define BCNROWSVARS 5    /**< Boundary condition's number of rows variable array position in vector. */

#define TIMER_NSECS 0    /**< To use nano second in Timer. */
#define TIMER_SECS 1    /**< To use second in Timer. */

#define G 9.81    /**< Gravitational acceleration. */
#define SQRTG 3.132091953    /**< Square root of Gravitational acceleration. */
#define EPS12 1e-12    /**< Tolerance e-12. */
#define FT3_TO_M3_FACTOR 0.028316847    /**< Factor to convert feet cube to meter cube. */
#define FT_TO_M_FACTOR 0.3048    /**< Factor to convert feet to meter. */
#define SEC_TO_HOUR_FACTOR 0.000277778    /**< Factor to convert second to hour. */
#define HOUR_TO_SEC_FACTOR 3600.0    /**< Factor to convert hour to second. */
#define MM_TO_M_FACTOR 0.001    /**< Factor to convert mili meter to meter. */

#define THREAD_BLOCK 256    /**< Thread block size to use in CUDA. */

#define TOTAL_TIME "total_time"    /**< Timer to get total runtime of the program. */
#define SIMULATION_TIME "simulation_time"    /**< Timer to get only the simulation time. */
#define COMPUTE_TIME "compute_time"    /**< Timer to get computation time. */
#define MPI_TIME "mpi_time"    /**< Timer to get all MPI operation time. */
#define IO_TIME "io_time"    /**< Timer to get time needed for outputting in file. */
#define RESIZE_TIME "resize_time"    /**< Timer to get time needed for resizing and re-balancing */
#define BALANCING_MPI_TIME "balancing_mpi_time"    /**< Timer to get time needed for resizing and re-balancing */


#define TYPE_STATIC "static"    /**< Domain decomposition type: static*/
#define TYPE_DYNAMIC "dynamic"    /**< Domain decomposition type: dynamic*/

#define RESET "\033[0m"     /**< Black Color */
#define RED "\033[31m"      /**< Red Color */
#define GREEN "\033[32m"    /**< Green Color */
#define YELLOW "\033[33m"   /**< Yellow Color */
#define BLUE "\033[34m"	    /**< Blue Color */
#define GRAY "\033[90m"	    /**< Gray Color */

#define OK GREEN << "[OK] " << RESET	    /**< Success Message */
#define WARN YELLOW << "[!!] " << RESET	    /**< Warning Message */
#define ERROR  RED << "[ERROR] " << RESET	/**< Error Message */
#define IN GRAY << "[..] " << RESET	        /**< Other Message 1 */
#define DASH BLUE << "[--] " << RESET	    /**< Other Message 2 */

#define WRITE_PERFORMANCE 0

#endif
