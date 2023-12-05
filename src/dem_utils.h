/** @file dem_utils.h
 *  @brief Header containing the DemFile class
 *
 *  This contains the subroutines and eventually any 
 *  macros, constants, etc. needed for DemFile class
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



#ifndef DEM_UTILS_H
#define DEM_UTILS_H

#include "matrix.h"

namespace DemFile
{
	template<class T>
	class dem_file : public Matrix::matrix<T>	/**< To process and store DEM data. It extends the base Matrix class. */
	{
	public:

/** @brief Constructor. Takes no argument.
*
*/
		dem_file<T>() : Matrix::matrix<T>() {}


/** @brief Constructor. Takes number of rows and columns as argument.
*
*  @param rows Number of rows
*  @param cols Number of columns
*/
		dem_file<T>(int rows, int cols) : Matrix::matrix<T>(rows, cols) {}
		
		
/** @brief Constructor. Takes a Matrix object as argument.
*
*  @param m Matrix object
*/		
		dem_file<T>(Matrix::matrix<T> const& m) : Matrix::matrix<T>(m) {}


/** @brief To get number of rows in DEM domain.
*
*  @return Number of rows
*/
		int get_nrows() const;
		
		
/** @brief To get number of columns in DEM domain.
*
*  @return Number of columns
*/		
		int get_ncols() const;
		
		
/** @brief To get the X coordinate of the origin
*
*  @return The X coordinate
*/
		T get_xll_corner() const;
		
		
/** @brief To get the Y coordinate of the origin
*
*  @return The Y coordinate
*/
		T get_yll_corner() const;
		
		
/** @brief To get the size of each cell
*
*  @return The cell size
*/		
		T get_cell_size() const;
		
		
/** @brief To get the default value if no data
*
*  @return No data value
*/
		int get_no_data_value() const;


/** @brief To set number of rows in DEM domain.
*
*  @param row Number of rows
*/
		void set_nrows(int row);
		
		
/** @brief To set number of columns in DEM domain.
*
*  @param col Number of columns
*/		
		void set_ncols(int col);
		
		
/** @brief To set X coordinate of the origin.
*
*  @param xll X coordinate
*/		
		void set_xll_corner(T xll);
		
		
/** @brief To set Y coordinate of the origin.
*
*  @param yll Y coordinate
*/		
		void set_yll_corner(T yll);
		
		
/** @brief To set size of a cell.
*
*  @param cell_size Cell size
*/		
		void set_cell_size(T cell_size);
		

/** @brief To set default value in case of no data.
*
*  @param no_data_value Deafult value
*/		
		void set_no_data_value(int no_data_value);


/** @brief Extract header information from a Ascii DEM file.
*
*  @param filename Ascii file name
*/
		void load_header_from_dem_file_ascii(std::string filename);
		
		
/** @brief Extract header information from a Binary DEM file.
*
*  @param filename Binary file name
*/		
		void load_header_from_dem_file_binary(std::string filename);

	private:
		int 
		nrows_,	/**< Number of rows in that domain or subdomain. */
		ncols_,	/**< Number of columns in that domain or subdomain. */
		no_data_value_;	/**< Default value in case of no data (not working now)*/
		
		T 
		xllcorner_, /**< X coordinate of the origin value. */
		yllcorner_,	/**< Y coordinate of the origin value. */
		cellsize_;	/**< Size of a cell. */
	};


	template<typename T>
	void dem_file<T>::load_header_from_dem_file_ascii(std::string filename)
	{
		std::ifstream ifs(filename.c_str());
		if (!ifs.good())
		{
			std::cerr << ERROR "Error reading file: " << filename << std::endl;
			exit(EXIT_FAILURE);

		}

		int line_num = 0;

		for (;;)
		{
			std::string line;
			std::getline(ifs, line);
			if (!ifs)
			break;

			line_num++;
			if (line_num > DEM_HEADER_SIZE)
			break;

			line = StringUtils::trim(line);
			std::vector<std::string> tokens = StringUtils::split(line, ' ');
			const char* value = (*(tokens.end() - 1)).c_str();

			if (!line.empty() && line[0] != '#')
			{
				switch (line_num)
				{
				case DEM_NCOLS_LINE:
				{
					this->ncols_ = atoi(value);
					break;
				}
				case DEM_NROWS_LINE:
				{
					this->nrows_ = atoi(value);
					break;
				}
				case DEM_XLL_CORNER_LINE:
				{
					this->xllcorner_ = atof(value);
					break;
				}
				case DEM_YLL_CORNER_LINE:
				{
					this->yllcorner_ = atof(value);
					break;
				}
				case DEM_CELL_SIZE_LINE:
				{
					this->cellsize_ = atof(value);
					break;
				}
				case DEM_NODATA_VALUE_LINE:
				{
					this->no_data_value_ = atoi(value);
					break;
				}
				default:
				{

				}
				}
			}
		}
		ifs.close();
	}

	
	template<typename T>
	void dem_file<T>::load_header_from_dem_file_binary(std::string filename)
	{
		std::ifstream ifs(filename.c_str(), ios::binary);
		if (!ifs.good())
		{
			std::cerr << ERROR "Error reading file: " << filename << std::endl;
			exit(EXIT_FAILURE);
		}

		T *arr = new T [DEM_HEADER_SIZE];

		ifs.read( (char*) arr, sizeof(T) * DEM_HEADER_SIZE );
		
		ifs.close();
		
		for(int i=0; i<DEM_HEADER_SIZE; i++)
		{
			int line_num = i+1;;
			T value = arr[i];
			
			switch (line_num)
			{
				case DEM_NCOLS_LINE:
				{
					this->ncols_ = (int)(value);
					break;
				}
				case DEM_NROWS_LINE:
				{
					this->nrows_ = (int)(value);
					break;
				}
				case DEM_XLL_CORNER_LINE:
				{
					this->xllcorner_ = value;
					break;
				}
				case DEM_YLL_CORNER_LINE:
				{
					this->yllcorner_ = value;
					break;
				}
				case DEM_CELL_SIZE_LINE:
				{
					this->cellsize_ = value;
					break;
				}
				case DEM_NODATA_VALUE_LINE:
				{
					this->no_data_value_ = (int)(value);
					break;
				}
				default:
				{

				}
			}	
		}
	}


	template<typename T>
	int dem_file<T>::get_nrows() const
	{
		return this->nrows_;
	}


	template<typename T>
	int dem_file<T>::get_ncols() const
	{
		return this->ncols_;
	}


	template<typename T>
	T dem_file<T>::get_cell_size() const
	{
		return this->cellsize_;
	}


	template<typename T>
	int dem_file<T>::get_no_data_value() const
	{
		return this->no_data_value_;
	}


	template<typename T>
	T dem_file<T>::get_xll_corner() const
	{
		return this->xllcorner_;
	}


	template<typename T>
	T dem_file<T>::get_yll_corner() const
	{
		return this->yllcorner_;
	}


	template<typename T>
	void dem_file<T>::set_nrows(int row)
	{
		this->nrows_ = row;
	}


	template<typename T>
	void dem_file<T>::set_ncols(int col)
	{
		this->ncols_ = col;
	}


	template<typename T>
	void dem_file<T>::set_xll_corner(T xll)
	{
		this->xllcorner_ = xll;
	}


	template<typename T>
	void dem_file<T>::set_yll_corner(T yll)
	{
		this->yllcorner_ = yll;
	}


	template<typename T>
	void dem_file<T>::set_cell_size(T cell_size)
	{
		this->cellsize_ = cell_size;
	}


	template<typename T>
	void dem_file<T>::set_no_data_value(int no_data_value)
	{
		this->no_data_value_ = no_data_value;
	}
}

#endif
