/** @file mpi_utils.h
 *  @brief Header containing the MpiUtils class
 *
 *  This contains the subroutines and eventually any 
 *  macros, constants, etc. needed for MpiUtils class
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



#ifndef MPI_UTILS_H
#define MPI_UTILS_H

#include "matrix.h"
#include "constants.h"

namespace MpiUtils
{
	
/** @brief It calculates each subdomain's number of rows and columns.
*
*  @param globalrows Row count of main domain
*  @param globalcols Column count of main domain
*  @param rank Current sub domain id
*  @param size Total number of sub domain
*  @return Row and column dimension
*/
	Constants::dims_t create_local_dims(int globalrows, int globalcols, int rank, int size);
	
	
	struct partition_data_t	/**< Structure to contain all subdomains row and column dimension. */
	{
		int size;	/**< Number of sub domains. */
		int rows;	/**< Number of rows in main domain. */
		int cols;	/**< Number of columns in main domain. */
		std::vector<Constants::dims_t > part_dims;	/**< Vector containing all subdomains dimension. */


/** @brief Constructor.
*
*/
		partition_data_t() : size(0), rows(0), cols(0) {}


/** @brief Constructor
*
*  @param s Number of sub domains
*  @param r Number of rows
*  @param c Number of columns
*/
		partition_data_t(int s, int r, int c) : size(s), rows(r), cols(c)
		{
			part_dims.assign(size, Constants::dims_t(0, 0));

			for (int p = 0; p < size; ++p)
			{
				part_dims[p] = create_local_dims(rows, cols, p, size);
			}
		}
	};
	
	
	Constants::dims_t create_local_dims(int globalrows, int globalcols, int rank, int size)
	{
		int halo = GHOST_CELL_PADDING;
		Constants::dims_t localgrid(halo * 2, halo * 2);

		localgrid.first += (globalrows / size);
		localgrid.second += globalcols;

		int rem = (globalrows % size);

		if (rank < rem && rem > 0)
		{
			localgrid.first++;
		}

		return localgrid;
	}


/** @brief It performs MPI halo exchanges between multiple MPI processes.
*
*  @param local Data to use for halo exchange
*  @param lrows Number of rows
*  @param lcols Number of columns
*  @param rank Current MPI process id
*  @param size Total number of MPI processes
*  @param type Data type (Only halo data/Full domain data)
*/
	template<typename T>
	void exchange(T* local, int lrows, int lcols, int rank, int size, int type);
	
	
/** @brief It performs initial domain scattering and partitioning between multiple MPI processes.
*
*  @param global Data of the whole domain
*  @param pd Partitioning information
*  @param rank Current MPI process id
*  @return Subdomain data
*/	
	template<typename T>
	Matrix::matrix<T> scatter_exchange(T* global, partition_data_t pd, int rank);
	
	
/** @brief It performs initial domain scattering and partitioning between multiple MPI processes for integer data type.
*
*  @param global Data of the whole domain
*  @param pd Partitioning information
*  @param rank Current MPI process id
*  @return Subdomain data
*/	
	Matrix::matrix<int> scatter_exchange_int(int* global, partition_data_t pd, int rank);


	template<typename T>
	void exchange(T* local, int lrows, int lcols, int rank, int size, int type)
	{
		int factor = 1;
		if (type == USE_HALO)
		{
			factor = lrows / 4;
		}

		int
		pr = (lcols * (lrows - 2 * factor)),
		lr = (lcols * (lrows - 1 * factor)),
		sr = (lcols * 1 * factor);

		int data_size = lcols * factor;

		MPI_Request send_request, recv_request;
		MPI_Status status;

		if (rank < size - 1)
		{
			MPI_Isend(&(local[pr]), data_size, MPI_DATA_TYPE, (rank + 1), 0, MPI_COMM_WORLD, &send_request);
		}

		if (rank > 0)
		{
			MPI_Irecv(&(local[0]), data_size, MPI_DATA_TYPE, (rank - 1), 0, MPI_COMM_WORLD, &recv_request);
		}

		if (rank < size - 1)
		{
			MPI_Wait(&send_request, &status);
		}

		if (rank > 0)
		{
			MPI_Wait(&recv_request, &status);
		}
		
		if (rank > 0)
		{
			MPI_Isend(&(local[sr]), data_size, MPI_DATA_TYPE, (rank - 1), 0, MPI_COMM_WORLD, &send_request);
		}

		if (rank < size - 1)
		{
			MPI_Irecv(&(local[lr]), data_size, MPI_DATA_TYPE, (rank + 1), 0, MPI_COMM_WORLD, &recv_request);
		}

		if (rank > 0)
		{
			MPI_Wait(&send_request, &status);
		}

		if (rank < size - 1)
		{
			MPI_Wait(&recv_request, &status);
		}


	}


	template<typename T>
	Matrix::matrix<T> scatter_exchange(T* global, partition_data_t pd, int rank)
	{
		int
		lrows = pd.part_dims[rank].first,
		lcols = pd.part_dims[rank].second,
		subsize = lrows * lcols;
		
		Matrix::matrix<T> sub_grid(lrows, lcols);
		sub_grid.zero_fill();
		MPI_Request send_request[pd.size-1], recv_request;

		if (rank == 0)
		{
			memcpy(sub_grid.get_address_at(0, 0), &global[0], subsize * sizeof(T));

			int row_pos = lrows - 2*GHOST_CELL_PADDING;
			
			for (int p = 1; p < pd.size; p++)
			{
				int
				lrows2 = pd.part_dims[p].first,
				lcols2 = pd.part_dims[p].second,
				subsize2 = lrows2 * lcols2;

				MPI_Isend(&global[row_pos*(long long)lcols2], subsize2, MPI_DATA_TYPE, p, 0, MPI_COMM_WORLD,&send_request[p-1]);
				row_pos += lrows2 - 2*GHOST_CELL_PADDING;
			}
		}
		else
		{
			MPI_Irecv(sub_grid.get_address_at(0, 0), subsize, MPI_DATA_TYPE, 0, 0, MPI_COMM_WORLD, &recv_request);
		}

		if (rank == 0)
		{
			MPI_Waitall(pd.size-1, send_request,  MPI_STATUS_IGNORE);
		}else{
			MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
		}

		return sub_grid;
	}


	Matrix::matrix<int> scatter_exchange_int(int* global, partition_data_t pd, int rank)
	{
		int
		lrows = pd.part_dims[rank].first,
		lcols = pd.part_dims[rank].second,
		subsize = lrows * lcols;

		Matrix::matrix<int> sub_grid(lrows, lcols);
		sub_grid.zero_fill_int();
		MPI_Request send_request[pd.size-1], recv_request;


		if (rank == 0)
		{
			memcpy(sub_grid.get_address_at(0, 0), &global[0], subsize * sizeof(int));
			
			int row_pos = lrows - 2*GHOST_CELL_PADDING;
			
			for (int p = 1; p < pd.size; p++)
			{
				int
				lrows2 = pd.part_dims[p].first,
				lcols2 = pd.part_dims[p].second,
				subsize2 = lrows2 * lcols2;

				MPI_Isend(&global[row_pos*(long long)lcols2], subsize2, MPI_INTEGER, p, 0, MPI_COMM_WORLD,&send_request[p-1]);
				row_pos += lrows2 - 2*GHOST_CELL_PADDING;
			}
		}
		else
		{
			MPI_Irecv(sub_grid.get_address_at(0, 0), subsize, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, &recv_request);
		}

		if (rank == 0)
		{
			MPI_Waitall(pd.size-1, send_request,  MPI_STATUS_IGNORE);
		}else{
			MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
		}

		return sub_grid;
	}
}

#endif
