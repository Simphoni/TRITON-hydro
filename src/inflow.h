/** @file inflow.h
 *  @brief Header containing the Hydrograph class
 *
 *  This contains the subroutines and eventually any 
 *  macros, constants, etc. needed for Hydrograph class
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



#ifndef INFLOW_H
#define INFLOW_H

#include "string_utils.h"

namespace Hydrograph
{
	template<class T>
	class hydrograph	/**< To process and store hydrograph input files. */
	{
	public:
		
/** @brief Constructor. Takes no argument.
*
*/	
		hydrograph();
		
		
/** @brief Constructor. Takes filename as an argument to construct the object.
*
*  @param filename Input file name
*/	
		hydrograph(std::string filename);


/** @brief It reads content from a hydrograph file and construct data.
*
*  @param filename Input file name
*/
		void load_from_file(std::string filename);


/** @brief To get all the contents in each rows of hydrograph file. 
*
*  @return All input rows.
*/
		std::vector<std::vector<T>> get_rows();
		
		
/** @brief It calculates flow value at a specific row index for a specific flow location number.
*
*  @param index Row index
*  @param source_num Flow location serial number
*  @return Flow value
*/		
		T get_flow_at(int index, int source_num);
		
		
/** @brief It calculates time at a specific row index.
*
*  @param index Row index
*  @return Time value
*/			
		T get_time_at(int index);


/** @brief To get number of inflow rows.
*
*  @return Inflow rows count
*/	
		int get_num_inflow_rows();
		
		
/** @brief To get number of inflows.
*
*  @return Inflows count
*/	
		int get_num_inflows();
		
		
/** @brief It converts all data time values from hour to second
*
*/			
		void convert_time_hr_to_secs();
		
		
/** @brief It converts all rate values from hour to second
*
*/	
		void convert_rate_hr_to_secs();
		
		
/** @brief It converts all mm values to m
*
*/			
		void convert_rate_mm_to_m();
		
		
/** @brief It sets number of flow rows.
*
*  @param rows Number of rows
*/			
		void set_num_flow_rows(int rows);
		
		
/** @brief It sets number of inflow locations.
*
*  @param sources inflow location count
*/			
		void set_num_sources(int sources);

	private:
		bool time_is_hours_;	/**< Flag to check time is in hours ot not. If true, time is in hours unit. */
		int flow_rows_;	  /**< Number of flow rows. */
		int num_sources_;	/**< Number of flow locations. */
		std::vector<std::vector<T> > data_;	  /**< Contains all the flow data. Inner vector contains data for each row. */
	};


	template<class T>
	hydrograph<T>::hydrograph()
	{
		time_is_hours_ = true;
		set_num_sources(1);
	}


	template<class T>
	hydrograph<T>::hydrograph(std::string filename)
	{
		time_is_hours_ = true;
		set_num_sources(1);
		load_from_file(filename);
		set_num_flow_rows(data_.size());
	}


	template<typename T>
	int hydrograph<T>::get_num_inflow_rows()
	{
		return flow_rows_;
	}


	template<typename T>
	int hydrograph<T>::get_num_inflows()
	{
		return num_sources_;
	}


	template<typename T>
	void hydrograph<T>::set_num_flow_rows(int rows)
	{
		flow_rows_ = rows;
	}


	template<typename T>
	void hydrograph<T>::set_num_sources(int sources)
	{
		num_sources_ = sources;
	}


	template<typename T>
	void hydrograph<T>::load_from_file(std::string filename)
	{
		std::ifstream ifs(filename.c_str());
		if (!ifs.good())
		{
			std::cerr << ERROR "Error reading file: " << filename << std::endl;
			exit(EXIT_FAILURE);
		}

		std::vector<std::vector<T> > hygdata;

		int line_num = 0;

		for (;;)
		{
			std::string line;
			std::getline(ifs, line);
			if (!ifs)
			break;

			line_num++;
			line = StringUtils::trim(line);

			std::vector<std::string> str_data = StringUtils::split(line, '%');

			if ((StringUtils::trim(str_data[0])).size() > 0)
			{
				std::vector<std::string> tokens = StringUtils::split(StringUtils::trim(str_data[0]), ',');
				std::vector<T> row = StringUtils::vecstr_to_vecflt<T>(tokens);

				data_.push_back(row);
			}
		}

		if(data_.size() > 0)
		{
			set_num_sources(data_[data_.size() - 1].size() - 1);
		}
		else
		{
			set_num_sources(0);
		}

	}


	template<typename T>
	std::vector<std::vector<T>> hydrograph<T>::get_rows()
	{
		return data_;
	}


	template<typename T>
	T hydrograph<T>::get_time_at(int index)
	{
		if (index >= static_cast<int>(data_.size()))
		{
			std::cerr << ERROR "Inflow or runoff index out of bounds: " << index << " out of " << data_.size() << std::endl;
			exit(EXIT_FAILURE);
		}

		return (data_.at(index)).at(0);
	}


	template<typename T>
	T hydrograph<T>::get_flow_at(int index, int source_num)
	{
		if (index >= static_cast<int>(data_.size()))
		{
			std::cerr << ERROR "Inflow or runoff index out of bounds: " << index << " out of " << data_.size() << std::endl;
			exit(EXIT_FAILURE);

		}

		return (data_.at(index)).at(source_num);
	}


	template<typename T>
	void hydrograph<T>::convert_time_hr_to_secs()
	{
		if (time_is_hours_)
		{
			typename std::vector<std::vector<T>>::iterator it = data_.begin();

			for (; it != data_.end(); it++)
			{
				(*it)[0] = it->at(0) * HOUR_TO_SEC_FACTOR;
			}

			time_is_hours_ = false;
		}
	}


	template<typename T>
	void hydrograph<T>::convert_rate_hr_to_secs()
	{
		for (int i = 0; i < get_num_inflows(); ++i)
		{
			for (int j = 0; j < get_num_inflow_rows(); j++)
			{
				data_[j][i + 1] /= HOUR_TO_SEC_FACTOR;
			}
		}

	}


	template<typename T>
	void hydrograph<T>::convert_rate_mm_to_m()
	{
		for (int i = 0; i < get_num_inflows(); ++i)
		{
			for (int j = 0; j < get_num_inflow_rows(); j++)
			{
				data_[j][i + 1] *= MM_TO_M_FACTOR;
			}
		}
	}
}

#endif
