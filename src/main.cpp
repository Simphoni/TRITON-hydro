/** @file main.cpp
 *  @brief Main file containing the driver
 *
 *  This contains the subroutines and eventually any 
 *  macros, constants, etc. needed for the driver
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



#include <iostream>
#include <string>
#include <float.h>
#include <utility>
#include <map>
#include <fstream>
#include <vector>
#include <sys/time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <dirent.h>
#include <cmath>
#include <cstring>
#include <omp.h>
#include "mpi.h"

using namespace std;

#include "constants.h"
#include "supertimer.h"
#include "string_utils.h"
#include "mpi_utils.h"
#include "config_utils.h"
#include "inflow.h"
#include "extbc.h"
#include "matrix.h"
#include "dem_utils.h"
#include "output.h"
#include "triton.h"
#include "kernels.h"


/** @brief Main function. This is the main function of the program.
*
*  @param argc Argument count
*  @param argv Pointer array which points to each argument passed to the program. The program runs with cfg filename and number of threads (only for OpenMP version)
*  @return 0
*/
int main(int argc, char* argv[])
{
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	Triton::triton<value_t> model(argc, argv);
	//initialize
	model.initialize(rank, size);
	//simulate
	model.simulate();

	MPI_Finalize();

	return 0;
}
