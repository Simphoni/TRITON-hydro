/** @file triton.h
 *  @brief Header containing the Triton class
 *
 *  This contains the subroutines and eventually any 
 *  macros, constants, etc. needed for Triton class
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



#ifndef TRITON_H
#define TRITON_H
#include "kernels.h"
#include "output.h"
#include "mpi_utils.h"

namespace Triton
{
	template<class T>
	class triton	/**< Main class to perform the simulation. */
	{
	public:
		
/** @brief Constructor.
*
*  @param argc Number of arguments
*  @param argv Arguments
*/	
		triton<T>(int argc, char* argv[]);
		
		
/** @brief Destructor. Releases any allocated memory.
*
*/
		~triton<T>();


/** @brief It initializes the simulation.
*
*  @param rank_ Subdomain id
*  @param size_ Number of subdomain
*/
		void initialize(int rank_, int size_);
		
		
/** @brief It starts the simulation. It is the main simulation fuction.
*
*/		
		void simulate();


	private:
		int rank;	/**< Current subdomain id */
		int size;	/**< Total sumber of subdomains */
		int rows;	/**< Number of rows in current subdomain */
		int cols;	/**< Number of columns in current subdomain */
		int org_rows;	/**< Number of rows in original domain without ghost cells */
		int org_cols;	/**< Number of columns in original domain without ghost cells */
		int num_of_src;	/**< Number of flow locations in current subdomain */
		int num_of_obs_points;	/**< Number of observations points in current subdomain */
		int num_of_extbc;	/**< Number of external boundary conditions */
		int num_extbc_cells;	/**< Number of cells in all external boundary conditions */
		int index_row_runoff;	/**< Index to keep track current runoff row id */
		int idx_low = 0;	/**< Index to keep track lower index id when updating flow locations */
		int checkpoint_id;	/**< Current checkpoint id */
		
		int host_src_pos_arr_size;	/**< Flow locations position container array size */
		int host_hyg_time_arr_size;	/**< Hydrograph time values container array size */
		int host_hyg_val_arr_size;	/**< Hydrograph flow values container array size */
		int host_bc_cells_size;	/**< Boundary condition cells container array size */
		int host_bc_vars_arr_size;	/**< Boundary condition variables container array size */
		int host_runoff_intensity_arr_size;	/**< Runoff intensity container array size */
		int host_reduce_dt_arr_sz;	/**< Time step size per cell container array size */
		int host_halo_arr_size;	/**< Halo cells container array size (water depth and discharges)*/
		int nbytes;	/**< Subdomain size in bytes */
		int nbytes_halo;	/**< Halo cells bundle size in bytes */
		
		std::vector<int> relative_obs_index; 	/**< Relative index position of observation cells per subdomain wrt to the global domain*/

		T simtime;	/**< Current simulation time */
		T cell_size;	/**< Cell size of the grid */
		T local_dt;	/**< Time step size in current subdomain */
		T global_dt;	/**< Time step size in the whole domain */
		T average_dt;	/**< Average time step size dince the last output */
		T init_dt;	/**< Time step size in case the domain is dry and the time step size is chosen dynamically*/
		
		std::string project_dir;	/**< Project directory */
		std::string cfg_content;	/**< Content of the input cfg file */
		std::string cfg_dir;	/**< Directory of input cfg file */

		SuperTimer::super_timer st;	/**< Time object to keep truct every custom timer */
		ConfigUtils::arguments<T> arglist;	/**< Object that holds all arguments and values from input cfg file */
		Hydrograph::hydrograph<T> hyg;	/**< Object that hold flow locations update data from hydrograph input file */
		Hydrograph::hydrograph<T> roff;	/**< Object that hold runoff input data */
		Constants::sources_list_t observation_cells;	/**< Local cell index information of all observation cells */
		Constants::sources_list_t observation_cells_global;	/**< Global cell index information of all observation cells */
		MpiUtils::partition_data_t pd;	/**< Partition information of all subdomains */
		MpiUtils::partition_data_t pd_aux;	/**< Auxuliar partition information of all subdomains in case of dynamic decomposition and parallel input */
		DemFile::dem_file<T> dem;	/**< Main domain's DEM file information and data */
		DemFile::dem_file<T> sub_dem;	/**< Current subdomain's DEM file information and data */
		Matrix::matrix<T> hin;	/**< Main domain's initial depth file data */
		Matrix::matrix<T> uin;	/**< Main domain's initial flux X data */
		Matrix::matrix<T> vin;	/**< Main domain's initial flux Y data */
		Matrix::matrix<T> nin;	/**< Main domain's initial manning data */
		Matrix::matrix<T> hot_hin;	/**< Main domain's water depth checkpoint data */
		Matrix::matrix<T> hot_qxin;	/**< Main domain's flux X checkpoint data */
		Matrix::matrix<T> hot_qyin;	/**< Main domain's flux Y checkpoint data */
		Matrix::matrix<T> sub_hin;	/**< Current subdomain's water depth data */
		Matrix::matrix<T> sub_qxin;	/**< Current subdomain's flux X data */
		Matrix::matrix<T> sub_qyin;	/**< Current subdomain's flux Y data */
		Matrix::matrix<T> sub_nin;	/**< Crrent subdomain's manning data */
		Matrix::matrix<T> sub_hot_hin;	/**< Current subdomain's water depth checkpoint data */
		Matrix::matrix<T> sub_hot_qxin;	/**< Current subdomain's discharge X checkpoint data */
		Matrix::matrix<T> sub_hot_qyin;	/**< Current subdomain's discharge Y checkpoint data */
		Matrix::matrix<T> sub_max_value_h;	/**< Current subdomain's max water depth of each cell */
		Matrix::matrix<T> hot_max_value_h;	/**< Main domain's max water depth checkpoint data */
		
		Matrix::matrix<int> rin;	/**< Main domain's runoff data */
		Matrix::matrix<int> sub_rin;	/**< Current subdomain's runoff data */
		
		int* host_src_pos_arr;	/**< Current subdomains flow locations position array */
		int* host_relative_bc_index;	/**< Boundary condition cells relative positions in current subdomain */
		int* host_bc_type;	/**< Boundary condition types of each boundary cells */
		int* host_bc_start_index;	/**< Separate boundary condition start indexes */
		int* host_bc_nrows_vars;	/**< Boundary condition variable of every rows */
		int* host_runoff_id_arr;	/**< Array that contains runoff ids */
		
		T* host_hyg_time_arr;	/**< Hydrograph time values container array */
		T* host_hyg_val_arr;	/**< Hydrograph flow values container array */
		T* host_extbc_var1_arr;	/**< Boundary condition variables container array */
		T* host_extbc_var2_arr;	/**< Boundary condition variables container array*/
		T* host_runoff_intensity_arr;	/**< Runoff intensity container array*/
		T* host_halo_arr;	/**< Halo cells container array*/
		T* host_sqrth_arr;	/**< Square root value array of every cells water depth */
		T* host_dt_values_arr;	/**< Time step size array of every cells */
		T* host_rhsh0;	/**< Intermediate raster array to hold partial water depth */
		T* host_rhsh1;	/**< Intermediate raster array to hold partial water depth */
		T* host_rhsqx0;	/**< Intermediate raster array to hold partial flux X */
		T* host_rhsqx1;	/**< Intermediate raster array to hold partial flux X */
		T* host_rhsqy0;	/**< Intermediate raster array to hold partial flux Y */
		T* host_rhsqy1;	/**< Intermediate raster array to hold partial flux Y */

		std::vector<T*> host_vec;	/**< Vector that contains all floating point array to use in simulation. */
		std::vector<int*> host_vec_int;	/**< Vector that contains all integer array to use in simulation. */

		Output::output<T> out; /**Object to manage output files. */ 

#ifdef ACTIVE_GPU
		cudaStream_t streams;	/**< Cuda stream */
		std::vector<T*> device_vec;	/**< Device vector that contains all floating point array to use in simulation. */
		std::vector<int*> device_vec_int;	/**< Device vector that contains all integer array to use in simulation. */
#endif

/** @brief This function is used to compute the time step size in the case the domain is dry and the dynamic CFL condition is set. It is computed according to the hydrograph interval, runoff interval and print interval divided by 100.
* 
*/
		void compute_init_dt();



/** @brief This function is used to calculate minimum time step size for each sub domain.
*
*/
		void compute_local_dt();
		
		
/** @brief This function is used to calculate minimum time step size between all sub domain.
*
*  @param print_id Current checkpoint id
*/		
		void compute_global_dt(int print_id);
		
		
/** @brief This function is used to compute next state. All main computation is done inside this function.
*
*/		
		void compute_new_state();
		
		
/** @brief It calculates a cell's column index.
*
*  @param src_x X location
*  @param xllc X coordinate of the origin
*  @param cell_size_ Cell size
*  @return The corresponding value
*/		
		int calc_src_col(T src_x, T xllc, T cell_size_);
		
		
/** @brief It calculates a cell's row index.
*
*  @param src_y Y location
*  @param yllc Y coordinate of the origin
*  @param cell_size_ Cell size
*  @param nrows Number of rows
*  @return The corresponding value
*/		
		int calc_src_row(T src_y, T yllc, T cell_size_, int nrows);
		
		
/** @brief This function reads from configuration file, process it for further use.
*
*  @param cfg_dir Cfg file directory
*  @param checkpoint_id Current checkpoint id
*/		
		void read_configuration(std::string cfg_dir, int checkpoint_id);
		
		
/** @brief This function reads all flow locations information.
*
*/			
		void read_inflows();


/** @brief This function reads the header of the dem files
*
*/		
		void read_header_dem_files();
		
		
/** @brief This function reads all matrix type data file in sequential mode.
*
*/		
		void read_matrix_files_sequential();

/** @brief This function reads all matrix type data file in parallel mode.
*
*/		
		void read_matrix_files_parallel();
		

/** @brief This function processes all flow locations and partition them in different subdomain.
*
*/		
		void process_source_locations();
		
		
/** @brief This function reads all observation cells information.
*
*/		
		void process_observation_cells();
		
		
/** @brief This function reads all boundary conditions and allow different subdomain capability.
*
*/
		void process_boundary_condition();
		
		
/** @brief This function processes all the matrix file and partition them into different sub domain.
*
*/		
		void partition_matrix_files();
		
		
/** @brief This function processes runoff if available.
*
*/		
		void process_runoff();
		
		
/** @brief This function creates additional arrays needed for simulation.
*
*/		
		void create_host_aux_vectors();
		
		
/** @brief This function creates host vector for simulation combining all needed arrays.
*
*/		
		void create_host_vectors();
		
		
/** @brief This function creates device vector for simulation combining all needed arrays.
*
*/		
		void create_device_vectors();


/** @brief This function creates new load balance domain decomposition.
*
*/		
		void new_domain_decomposition();

		/** @brief This function gathers all the information from dynamic partition, creates the new partition based on previous MPI times and broadcasts the information from process 0 to the rest of the processes.
*
*/		
		int MPI_time_based_domain_decomposition();

		/** @brief This function resets previous arrays.
*
*/		
		void reset_arrays();


/** @brief This function processes all the matrix file and partition them into different sub domain for dynamic balancing
*
*/		
		void partition_matrix_files_dynamic();



	};


	template<class T>
	triton<T>::triton(int argc, char* argv[])
	{
		st = SuperTimer::super_timer();
		st.start(TOTAL_TIME);
		project_dir = ConfigUtils::get_root_dir(argv[0]);
		if(argc > 1)
		{
			cfg_dir = std::string(argv[1]);
		}
		
#ifdef ACTIVE_OMP
		int threads = 1;
		if(argc > 2)
		{
			threads = atoi(argv[2]);
		}
		omp_set_num_threads(threads);
#endif
		
		checkpoint_id = 0;
		if(argc > 3)
		{
			checkpoint_id = atoi(argv[3]);
		}

	}


	template<typename T>
	void triton<T>::initialize(int rank_, int size_)
	{
		rank = rank_;
		size = size_;

		read_configuration(cfg_dir, checkpoint_id);
		simtime = arglist.sim_start_time;
		
		read_inflows();

		read_header_dem_files();
		
		if(strcmp(arglist.input_option.c_str(), "SEQ")==0){
			read_matrix_files_sequential();
		}

		//the first partitioning at the beginning of the simulation is homogeneous (static) unless there is checkpoint, more than 1 processes and dynamic partition
		pd = MpiUtils::partition_data_t(size, org_rows, org_cols);
		if (arglist.checkpoint_id > 0 && size > 1 && strcmp(arglist.domain_decomposition.c_str(), TYPE_DYNAMIC)==0){
			int *dyn_rows = new int[size];
			if(rank==0){
				ConfigUtils::read_and_parse_checkpoint_partition(project_dir, dyn_rows, arglist.checkpoint_id);
			}
			MPI_Bcast(dyn_rows, size, MPI_INT, 0, MPI_COMM_WORLD); 
			for(int i=0;i<pd.size;i++){
				pd.part_dims[i].first=dyn_rows[i]+2*GHOST_CELL_PADDING;
			}
		}

		process_source_locations();
		process_observation_cells();
		
		process_boundary_condition();

		if(strcmp(arglist.input_option.c_str(), "SEQ")==0){
			partition_matrix_files();
		}else{
			//we need pd_aux to read the dem, the rmap and the mann files with a "regular" decomposition and then convert them to the dynamic 
			pd_aux = MpiUtils::partition_data_t(size, org_rows, org_cols);
			read_matrix_files_parallel();
		}	

		process_runoff();

		create_host_aux_vectors();
		
		create_host_vectors();
		create_device_vectors();
	}
	
	
	template<typename T>
	void triton<T>::read_configuration(std::string cfg_dir, int checkpoint_id)
	{
		if (rank == 0){
			std::cerr << IN "Reading configuration file" << std::endl;
		}

		std::string cfg_path = project_dir + "/" + INPUT_DIR + "/" + CFG_DIR + "/" + DEFAULT_CFG;
		
		if(!(cfg_dir.empty()) && !(StringUtils::is_numeric(cfg_dir)))
		{
			cfg_path = cfg_dir;
		}

		if (checkpoint_id > 0)
		{
			cfg_path = project_dir + "/" + OUTPUT_DIR + "/" + CFG_DIR + "/config_" + to_string(checkpoint_id) + ".cfg";
		}

		cfg_content = ConfigUtils::file_content_to_string(cfg_path);
		arglist = ConfigUtils::get_args<T>(cfg_content);

		if (rank == 0){
			std::cerr << OK "Configuration file read" << std::endl;
			if (arglist.num_sources == 0)
			{
				std::cerr << DASH "No sources defined" << std::endl;
			}else{
				std::cerr << DASH << arglist.num_sources << " sources defined" << std::endl;
			}

			if (arglist.num_runoffs == 0)
			{
				std::cerr << DASH "No runoff defined" << std::endl;
			}else{
				std::cerr << DASH << arglist.num_runoffs << " runoffs defined" << std::endl;
			}

			if (arglist.num_extbc == 0)
			{
				std::cerr << DASH "No external boundary conditions defined" << std::endl;
			}else{
				std::cerr << DASH << arglist.num_extbc << " external boundary conditions defined" << std::endl;
			}

			if (arglist.time_series_flag == 0)
			{
				std::cerr << DASH "No observation points defined" << std::endl;
			}else{
				std::cerr << DASH << arglist.observation_x_loc.size() << " observation points defined" << std::endl;
			}
		}
	}
	
	
	template<typename T>
	void triton<T>::read_inflows()
	{
		if(arglist.num_sources > 0)
		{		
			if (rank == 0){
				std::cerr << IN "Reading and processing source hydrographs" << std::endl;
			}
			hyg = Hydrograph::hydrograph<T>(arglist.hydrograph_filename);
			hyg.convert_time_hr_to_secs();
			if (rank == 0){
				std::cerr << OK "Sources set" << std::endl;
			}
		}
		
		if (arglist.num_runoffs > 0)
		{
			if (rank == 0){
				std::cerr << IN "Reading and processing runoff rates" << std::endl;
			}
			roff = Hydrograph::hydrograph<T>(arglist.runoff_filename);
			roff.convert_time_hr_to_secs();
			roff.convert_rate_hr_to_secs();
			roff.convert_rate_mm_to_m();
			if (rank == 0){
				std::cerr << OK "Runoff rates set" << std::endl;
			}
		}
	}
	
	
	template<typename T>
	void triton<T>::read_header_dem_files()
	{
		if(strcmp(arglist.input_option.c_str(), "SEQ")==0){

			if (strcmp(arglist.input_format.c_str(), "ASC") == 0)
			{
				dem.load_header_from_dem_file_ascii(arglist.dem_filename);
			}
			else
			{
				dem.load_header_from_dem_file_binary(arglist.dem_filename);
			}
		}else{
			//same subroutine to read the header, but from the header file instead of the full dem file
			//the header will be read from ASCII always, regardless of the input_format
			dem.load_header_from_dem_file_ascii(arglist.header_filename);
		}
		
		org_rows = dem.get_nrows();
		org_cols = dem.get_ncols();
		cell_size = dem.get_cell_size();
	}

	template<typename T>
	void triton<T>::read_matrix_files_sequential()
	{

		if (rank == 0)
		{
			if (strcmp(arglist.input_format.c_str(), "ASC") == 0)
			{
				dem.load_from_ascii_file(org_rows, org_cols, arglist.dem_filename, DEM_HEADER_SIZE);
			}
			else
			{
				dem.load_from_binary_file(org_rows, org_cols, arglist.dem_filename, DEM_HEADER_SIZE);
			}
			dem.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
			dem.set_nrows(dem.get_num_rows());
			dem.set_ncols(dem.get_num_cols());
			
			if(!arglist.open_boundaries){
				dem.set_infinite_walls();
			}else{
				dem.copy_value_into_ghost_cells();
			}

			if(!arglist.n_infile.empty())
			{
				if (strcmp(arglist.input_format.c_str(), "ASC") == 0)
				{
					nin.load_from_ascii_file(org_rows, org_cols, arglist.n_infile);
				}
				else
				{
					nin.load_from_binary_file(org_rows, org_cols, arglist.n_infile);
				}
			}
			else
			{
				nin.resize(org_rows, org_cols);
				nin.zero_fill();
				nin += arglist.const_mann;
			}
			nin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
			nin.copy_value_into_ghost_cells();
			nin.square();
			
			if (arglist.h_infile.size() > 0)
			{
				if (strcmp(arglist.input_format.c_str(), "ASC") == 0)
				{
					hin.load_from_ascii_file(org_rows, org_cols, arglist.h_infile);
				}
				else
				{
					hin.load_from_binary_file(org_rows, org_cols, arglist.h_infile);
				}
				hin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
				
				if(arglist.open_boundaries){				
					hin.copy_value_into_ghost_cells();
				}
			}
			if (arglist.qx_infile.size() > 0)
			{
				if (strcmp(arglist.input_format.c_str(), "ASC") == 0)
				{
					uin.load_from_ascii_file(org_rows, org_cols, arglist.qx_infile);
				}
				else
				{
					uin.load_from_binary_file(org_rows, org_cols, arglist.qx_infile);
				}
				uin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
				
				if(arglist.open_boundaries){				
					uin.copy_value_into_ghost_cells();
				}

			}
			if (arglist.qy_infile.size() > 0)
			{
				if (strcmp(arglist.input_format.c_str(), "ASC") == 0)
				{
					vin.load_from_ascii_file(org_rows, org_cols, arglist.qy_infile);
				}
				else
				{
					vin.load_from_binary_file(org_rows, org_cols, arglist.qy_infile);
				}
				vin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
				
				if(arglist.open_boundaries){				
					vin.copy_value_into_ghost_cells();
				}

			}

			if (arglist.runoff_map.size() > 0)
			{
				if (strcmp(arglist.input_format.c_str(), "ASC") == 0)
				{
					rin.load_from_ascii_file(org_rows, org_cols, arglist.runoff_map);
				}
				else
				{
					rin.load_from_binary_file(org_rows, org_cols, arglist.runoff_map);
				}
				rin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0);
			}
			
			if (arglist.checkpoint_id > 0 && strcmp(arglist.output_option.c_str(), "SEQ") == 0)
			{
				if (rank == 0){
					std::cerr << IN "Reading checkpoint files" << std::endl;
				}
				string temp_num(to_string(arglist.checkpoint_id));
				if (arglist.checkpoint_id < 10)
				{
					temp_num = "0" + temp_num;
				}

				string filedirH(project_dir + "/" + OUTPUT_DIR + "/" + BIN_DIR + "/H_" + temp_num + "_00.out");
				hot_hin.load_from_binary_file(org_rows, org_cols, filedirH);
				hot_hin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
				if(arglist.open_boundaries){				
					hot_hin.copy_value_into_ghost_cells();
				}

				string filedirQX(project_dir + "/" + OUTPUT_DIR + "/" + BIN_DIR + "/QX_" + temp_num + "_00.out");
				hot_qxin.load_from_binary_file(org_rows, org_cols, filedirQX);
				hot_qxin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
				if(arglist.open_boundaries){				
					hot_qxin.copy_value_into_ghost_cells();
				}


				string filedirQY(project_dir + "/" + OUTPUT_DIR + "/" + BIN_DIR + "/QY_" + temp_num + "_00.out");
				hot_qyin.load_from_binary_file(org_rows, org_cols, filedirQY);
				hot_qyin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
				
				if(arglist.open_boundaries){				
					hot_qyin.copy_value_into_ghost_cells();
				}

				
				if (arglist.max_value_print_option.size() > 0)	
				{
					string filedirMaxH(project_dir + "/" + OUTPUT_DIR + "/" + BIN_DIR + "/MH_" + temp_num + "_00.out");
					hot_max_value_h.load_from_binary_file(org_rows, org_cols, filedirMaxH);
					hot_max_value_h.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
					if(arglist.open_boundaries){				
						hot_max_value_h.copy_value_into_ghost_cells();
					}

				}
				
				if (rank == 0){
					std::cerr << OK "Checkpoint files read" << std::endl;
				}

			}
		}
	}


	template<typename T>
	void triton<T>::read_matrix_files_parallel()
	{

		int lrows = pd.part_dims[rank].first - 2 * GHOST_CELL_PADDING;
		int lcols = pd.part_dims[rank].second - 2 * GHOST_CELL_PADDING;

		//regular decomposition, just from dem, mann and rmap. In the usual case they match lrows and lcols. They will only differ for dynamic domain decomposition
		int lrows1 = pd_aux.part_dims[rank].first - 2 * GHOST_CELL_PADDING;
		int lcols1 = pd_aux.part_dims[rank].second - 2 * GHOST_CELL_PADDING;
		
		string temp_rank(to_string(rank));
		if (rank < 10)
		{
			temp_rank = "0" + temp_rank;
		}

		string filedir_dem(arglist.dem_filename + "_" + temp_rank + ".dem");
		if (strcmp(arglist.input_format.c_str(), "ASC") == 0)
		{
			sub_dem.load_from_ascii_file(lrows1, lcols1, filedir_dem, 0); //0 is because there is no header_dem_size in parallel reading. The header is in a separate file		
		}
		else
		{
			sub_dem.load_from_binary_file(lrows1, lcols1, filedir_dem); //the binary reading always have two values at the beginnig with the number of rows and columns.
		}

		sub_dem.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
		
		if (arglist.checkpoint_id > 0 && size > 1 && strcmp(arglist.domain_decomposition.c_str(), TYPE_DYNAMIC)==0){
			//we need to call out.init to have all the information in the struct out that is used afterwards. we need (lrows1,lcols1) plus the ghost cells as arguments 
			//since the subroutine workis with the full subdomain (including ghost cells)
			out.init(lrows1+2 * GHOST_CELL_PADDING, lcols1+2 * GHOST_CELL_PADDING, rank, size, project_dir, arglist.outfile_pattern, arglist.time_series_flag, cfg_content, arglist.output_option);
			if (rank == 0)
			{
				MPI_Gatherv(sub_dem.get_address_at(0, 0), out.cur_proc_data_size, MPI_DATA_TYPE, out.total_data_arr, out.recvcounts, out.displs, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);
			}
			else
			{
				MPI_Gatherv(sub_dem.get_address_at(GHOST_CELL_PADDING, 0), out.cur_proc_data_size, MPI_DATA_TYPE, out.total_data_arr, out.recvcounts, out.displs, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);
			}
			//now in out.total_data_arr we have the full dem that we have to partition according to the last state in the dynamic decomposition
			sub_dem.resize(1,1);
			sub_dem = MpiUtils::scatter_exchange(out.total_data_arr, pd, rank);
		}
		
		rows = sub_dem.get_num_rows();
		cols = sub_dem.get_num_cols();
		
		sub_dem.set_nrows(sub_dem.get_num_rows());
		sub_dem.set_ncols(sub_dem.get_num_cols());
		sub_dem.set_cell_size(dem.get_cell_size());
		sub_dem.set_xll_corner(dem.get_xll_corner());
		sub_dem.set_yll_corner(dem.get_yll_corner());
		sub_dem.set_no_data_value(dem.get_no_data_value());

		if(!arglist.open_boundaries){
			sub_dem.set_infinite_walls();
		}else{
			sub_dem.copy_value_into_ghost_cells();
		}

		MpiUtils::exchange(sub_dem.begin(), rows, cols, rank, size, USE_MATRIX);
		MPI_Barrier(MPI_COMM_WORLD);
		
		//modify ghost values that contains any sort of external boundary condition
		if (num_of_extbc > 0 && num_extbc_cells > 0){
			for(int i=0;i<num_extbc_cells;i++){
				int ii = host_relative_bc_index[i];
				int ix = (ii / cols);	//row id
				int iy = (ii % cols);	//col id
				bool
				is_top = (ix == 1),
				is_btm = (ix == rows - 2),
				is_lt = (iy == 1),
				is_rt = (iy == cols - 2);

				if(is_lt){ //west
					sub_dem.set_value(ii-1, sub_dem.get_value(ii));
				}
				if(is_rt){ //east
					sub_dem.set_value(ii+1, sub_dem.get_value(ii));
				}
				if (rank == 0 && is_top){ //north
					sub_dem.set_value(ii-cols, sub_dem.get_value(ii));
				}
				if (rank == size - 1 && is_btm) //south
				{
					sub_dem.set_value(ii+cols, sub_dem.get_value(ii));
				}
			
			}
			
		}


		if(!arglist.n_infile.empty())
		{
			string filedir_nin(arglist.n_infile + "_" + temp_rank + ".mann");
			
			if (strcmp(arglist.input_format.c_str(), "ASC") == 0)
			{
				sub_nin.load_from_ascii_file(lrows1, lcols1, filedir_nin);
			}
			else
			{
				sub_nin.load_from_binary_file(lrows1, lcols1, filedir_nin);
			}
		}
		else
		{
			//if there is no manning file we still have to build the sub_nin matrices with lrows1,lcols1 . We can just do it with lrows,lcols and do nothing for the dynamic decomposition 
			sub_nin.resize(lrows, lcols);
			sub_nin.zero_fill();
			sub_nin += arglist.const_mann;
		}
		sub_nin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);

		if (arglist.checkpoint_id > 0 && size > 1 && strcmp(arglist.domain_decomposition.c_str(), TYPE_DYNAMIC)==0 && !arglist.n_infile.empty()){
			//note that we do this just in the case of the user provided with a mann file. If it is a constant number we don't have to do it since we assigned correctly the size using lrows,lcols
			if (rank == 0)
			{
				MPI_Gatherv(sub_nin.get_address_at(0, 0), out.cur_proc_data_size, MPI_DATA_TYPE, out.total_data_arr, out.recvcounts, out.displs, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);
			}
			else
			{
				MPI_Gatherv(sub_nin.get_address_at(GHOST_CELL_PADDING, 0), out.cur_proc_data_size, MPI_DATA_TYPE, out.total_data_arr, out.recvcounts, out.displs, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);
			}
			//now in out.total_data_arr we have the full nin that we have to partition according to the last state in the dynamic decomposition
			sub_nin.resize(1,1);
			sub_nin = MpiUtils::scatter_exchange(out.total_data_arr, pd, rank);

		}

		sub_nin.copy_value_into_ghost_cells();
		sub_nin.square();

		MpiUtils::exchange(sub_nin.begin(), rows, cols, rank, size, USE_MATRIX);
		MPI_Barrier(MPI_COMM_WORLD);


		if (arglist.h_infile.size() > 0)
		{
			string filedir_hin(arglist.h_infile + "_" + temp_rank + ".out");

			if (strcmp(arglist.input_format.c_str(), "ASC") == 0)
			{
				sub_hin.load_from_ascii_file(lrows, lcols, filedir_hin);
			}
			else
			{
				sub_hin.load_from_binary_file(lrows, lcols, filedir_hin);
			}

			sub_hin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
			
			if(arglist.open_boundaries){				
				sub_hin.copy_value_into_ghost_cells();
			}
			
			MpiUtils::exchange(sub_hin.begin(), rows, cols, rank, size, USE_MATRIX);
			MPI_Barrier(MPI_COMM_WORLD);
		
		}

		if (arglist.qx_infile.size() > 0)
		{
			string filedir_qxin(arglist.qx_infile + "_" + temp_rank + ".out");

			if (strcmp(arglist.input_format.c_str(), "ASC") == 0)
			{
				sub_qxin.load_from_ascii_file(lrows, lcols, filedir_qxin);
			}
			else
			{
				sub_qxin.load_from_binary_file(lrows, lcols, filedir_qxin);
			}
			sub_qxin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
			
			if(arglist.open_boundaries){				
				sub_qxin.copy_value_into_ghost_cells();
			}
			MpiUtils::exchange(sub_qxin.begin(), rows, cols, rank, size, USE_MATRIX);
			MPI_Barrier(MPI_COMM_WORLD);


		}

		if (arglist.qy_infile.size() > 0)
		{
			string filedir_qyin(arglist.qy_infile + "_" + temp_rank + ".out");

			if (strcmp(arglist.input_format.c_str(), "ASC") == 0)
			{
				sub_qyin.load_from_ascii_file(lrows, lcols, filedir_qyin);
			}
			else
			{
				sub_qyin.load_from_binary_file(lrows, lcols, filedir_qyin);
			}
			sub_qyin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
			
			if(arglist.open_boundaries){				
				sub_qyin.copy_value_into_ghost_cells();
			}
			MpiUtils::exchange(sub_qyin.begin(), rows, cols, rank, size, USE_MATRIX);
			MPI_Barrier(MPI_COMM_WORLD);

		}

		if (arglist.runoff_map.size() > 0)
		{
			string filedir_rmap(arglist.runoff_map + "_" + temp_rank + ".rmap");

			if (strcmp(arglist.input_format.c_str(), "ASC") == 0)
			{
				sub_rin.load_from_ascii_file(lrows1, lcols1, filedir_rmap);
			}
			else
			{
				sub_rin.load_from_binary_file(lrows1, lcols1, filedir_rmap);
			}
			sub_rin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0);

			if (arglist.checkpoint_id > 0 && size > 1 && strcmp(arglist.domain_decomposition.c_str(), TYPE_DYNAMIC)==0){
				
				if (rank == 0)
				{
					MPI_Gatherv(sub_rin.get_address_at(0, 0), out.cur_proc_data_size, MPI_INTEGER, out.total_data_arr_int, out.recvcounts, out.displs, MPI_INTEGER, 0, MPI_COMM_WORLD);
				}
				else
				{
					MPI_Gatherv(sub_rin.get_address_at(GHOST_CELL_PADDING, 0), out.cur_proc_data_size, MPI_INTEGER, out.total_data_arr_int, out.recvcounts, out.displs, MPI_INTEGER, 0, MPI_COMM_WORLD);
				}
				sub_rin.resize(1,1);
				sub_rin = MpiUtils::scatter_exchange_int(out.total_data_arr_int, pd, rank);

			}

			sub_rin.copy_value_into_ghost_cells();

			//not neccesary to exchange. Otherwise a new function should be done because MpiUtils::exchange works with doubles
			//MpiUtils::exchange(sub_rin.begin(), rows, cols, rank, size, USE_MATRIX);
			//MPI_Barrier(MPI_COMM_WORLD);
		}
		

		if (arglist.checkpoint_id > 0)
		{
			if (rank == 0){
				std::cerr << IN "Reading checkpoint files" << std::endl;
			}
			string temp_num(to_string(arglist.checkpoint_id));
			if (arglist.checkpoint_id < 10)
			{
				temp_num = "0" + temp_num;
			}

			string filedirH(project_dir + "/" + OUTPUT_DIR + "/" + BIN_DIR + "/H_" + temp_num + "_" + temp_rank + ".out");
			sub_hot_hin.load_from_binary_file(lrows, lcols, filedirH);
			sub_hot_hin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
			if(arglist.open_boundaries){				
				sub_hot_hin.copy_value_into_ghost_cells();
			}

			string filedirQX(project_dir + "/" + OUTPUT_DIR + "/" + BIN_DIR + "/QX_" + temp_num + "_" + temp_rank + ".out");
			sub_hot_qxin.load_from_binary_file(lrows, lcols, filedirQX);
			sub_hot_qxin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
			if(arglist.open_boundaries){				
				sub_hot_qxin.copy_value_into_ghost_cells();
			}


			string filedirQY(project_dir + "/" + OUTPUT_DIR + "/" + BIN_DIR + "/QY_" + temp_num + "_" + temp_rank + ".out");
			sub_hot_qyin.load_from_binary_file(lrows, lcols, filedirQY);
			sub_hot_qyin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
			
			if(arglist.open_boundaries){				
				sub_hot_qyin.copy_value_into_ghost_cells();
			}

			MpiUtils::exchange(sub_hot_hin.begin(), rows, cols, rank, size, USE_MATRIX);
			MPI_Barrier(MPI_COMM_WORLD);
			MpiUtils::exchange(sub_hot_qxin.begin(), rows, cols, rank, size, USE_MATRIX);
			MPI_Barrier(MPI_COMM_WORLD);
			MpiUtils::exchange(sub_hot_qyin.begin(), rows, cols, rank, size, USE_MATRIX);
			MPI_Barrier(MPI_COMM_WORLD);

			
			if (arglist.max_value_print_option.size() > 0)	
			{
				string filedirMaxH(project_dir + "/" + OUTPUT_DIR + "/" + BIN_DIR + "/MH_" + temp_num + "_" + temp_rank + ".out");
				sub_max_value_h.load_from_binary_file(lrows, lcols, filedirMaxH);
				sub_max_value_h.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
				if(arglist.open_boundaries){				
					sub_max_value_h.copy_value_into_ghost_cells();
				}

			}
			
			if (rank == 0){
				std::cerr << OK "Checkpoint files read" << std::endl;
			}

		}

		//overwrite hin, qxin, qyin 
		if (arglist.checkpoint_id > 0)
		{
			sub_hin = sub_hot_hin;
			sub_qxin = sub_hot_qxin;
			sub_qyin = sub_hot_qyin;
			
			if (arglist.max_value_print_option.size() <= 0)	
			{
				sub_max_value_h = sub_hin;
			}
		}
		else
		{
			Matrix::matrix<T> values(rows, cols);
			values.zero_fill();

			if (arglist.h_infile.size() <= 0)
			{
				sub_hin = values;
			}
			if (arglist.qx_infile.size() <= 0)
			{
				sub_qxin = values;
			}
			if (arglist.qy_infile.size() <= 0)
			{
				sub_qyin = values;
			}
			
			sub_max_value_h = sub_hin;
		}

		if (rank == 0){
			std::cerr << OK "Data read in parallel" << std::endl;
		}


	}



	
	
	template<typename T>
	void triton<T>::process_source_locations()
	{
		num_of_src = 0;
		std::vector<int> src_rows, src_cols;
		Constants::sources_list_t source_cells;
		std::vector<int> relative_src_index;
		
		if(arglist.num_sources > 0)
		{
			std::vector<T> src_x = arglist.src_x_loc;
			std::vector<T> src_y = arglist.src_y_loc;

			int num_sources = arglist.num_sources;
			src_rows.assign(num_sources, 0);
			src_cols.assign(num_sources, 0);

			for (int i = 0; i < num_sources; ++i)
			{
				src_cols[i] = calc_src_col(src_x[i], dem.get_xll_corner(), dem.get_cell_size());
				src_rows[i] = calc_src_row(src_y[i], dem.get_yll_corner(), dem.get_cell_size(), org_rows);
				
				if(src_cols[i] >= org_cols || src_rows[i] >= org_rows || src_cols[i]<0 || src_rows[i]<0){
					std::cerr << ERROR "Source " << i+1  << " is out of bounds" << std::endl;
					exit(EXIT_FAILURE);
				}
			}

			for (int i = 0; i < arglist.num_sources; ++i)
			{
				int srank = 0;
				int prev_rows_sum = 0;

				if(size > 1){
					int source_row = src_rows[i];
					int rows_sum = pd.part_dims[0].first - 2 * GHOST_CELL_PADDING;
					
					if(source_row >= rows_sum){
						for(int j=1; j<size; j++){
							prev_rows_sum = rows_sum;
							rows_sum += pd.part_dims[j].first - 2 * GHOST_CELL_PADDING;
							if(source_row < rows_sum){
								srank = j;
								break;
							}
						}
					}
				}
				src_rows[i] = src_rows[i] - prev_rows_sum + GHOST_CELL_PADDING;
				src_cols[i] = src_cols[i] + GHOST_CELL_PADDING;
				
				if (rank == srank)
				{
					relative_src_index.push_back(i);
					
					std::pair<int, int> scell(src_rows[i], src_cols[i]);
					source_cells.push_back(scell);

					num_of_src++;
				}
			}
		}
		
		host_src_pos_arr_size = num_of_src;
		if (host_src_pos_arr_size == 0)
		{
			host_src_pos_arr_size = 1;
		}
		host_src_pos_arr = new int[host_src_pos_arr_size];
		for (int i = 0; i < host_src_pos_arr_size; i++)
		{
			host_src_pos_arr[i] = 0;
		}

		for (int i = 0; i<num_of_src; i++)
		{
			std::pair<int, int> pair = source_cells[i];
			host_src_pos_arr[i] = pair.first*(org_cols+2*GHOST_CELL_PADDING) + pair.second;
		}

		host_hyg_time_arr_size = 1;
		if(num_of_src > 0)
		{
			host_hyg_time_arr_size = hyg.get_num_inflow_rows();
		}
		
		host_hyg_time_arr = new T[host_hyg_time_arr_size];
		
		host_hyg_val_arr_size = host_hyg_time_arr_size*host_src_pos_arr_size;
		host_hyg_val_arr = new T[host_hyg_val_arr_size];
		if(num_of_src > 0)
		{
			for (int i = 0; i < host_hyg_time_arr_size; i++)
			{
				host_hyg_time_arr[i] = hyg.get_time_at(i);
				
				for (int j = 0; j < num_of_src; j++)
				{
					int src_number = relative_src_index[j] + 1;
					host_hyg_val_arr[i*num_of_src + j] = hyg.get_flow_at(i, src_number);
				}
			}
		}
		else
		{
			host_hyg_time_arr[0] = 0.0;
			host_hyg_val_arr[0] = 0.0;
		}
	}
	
	
	template<typename T>
	void triton<T>::process_observation_cells()
	{

		num_of_obs_points = 0;

		if (arglist.time_series_flag)
		{
			std::vector<T> observation_x = arglist.observation_x_loc;
			std::vector<T> observation_y = arglist.observation_y_loc;

			int num_observation_loc_global = observation_x.size();
			std::vector<int> observation_rows, observation_cols;

			observation_rows.assign(num_observation_loc_global, 0);
			observation_cols.assign(num_observation_loc_global, 0);

			for (int i = 0; i < num_observation_loc_global; ++i)
			{
				observation_cols[i] = calc_src_col(observation_x[i], dem.get_xll_corner(), dem.get_cell_size());
				observation_rows[i] = calc_src_row(observation_y[i], dem.get_yll_corner(), dem.get_cell_size(), org_rows);
				
				if(observation_cols[i] >= org_cols || observation_rows[i] >= org_rows || observation_cols[i]<0 || observation_rows[i]<0){
					std::cerr << ERROR "Observation " << i+1  << " is out of bounds" << std::endl;
				}
				std::pair<int, int> scell(observation_rows[i]+ GHOST_CELL_PADDING, observation_cols[i]+ GHOST_CELL_PADDING);
				observation_cells_global.push_back(scell);

			}


			for (int i = 0; i < num_observation_loc_global; ++i)
			{
				int srank = 0;
				int prev_rows_sum = 0;

				if(size > 1){
					int observation_row = observation_rows[i];
					int rows_sum = pd.part_dims[0].first - 2 * GHOST_CELL_PADDING;
					
					if(observation_row >= rows_sum){
						for(int j=1; j<size; j++){
							prev_rows_sum = rows_sum;
							rows_sum += pd.part_dims[j].first - 2 * GHOST_CELL_PADDING;
							if(observation_row < rows_sum){
								srank = j;
								break;
							}
						}
					}
				}
				observation_rows[i] = observation_rows[i] - prev_rows_sum + GHOST_CELL_PADDING;
				observation_cols[i] = observation_cols[i] + GHOST_CELL_PADDING;
				
				if (rank == srank)
				{
					relative_obs_index.push_back(i);
					std::pair<int, int> scell(observation_rows[i], observation_cols[i]);
					observation_cells.push_back(scell);
					num_of_obs_points++;
				}
			}

		}

	}
	
	
	template<typename T>
	void triton<T>::process_boundary_condition()
	{
		int num_extbc = 1;
		if(arglist.num_extbc > 0) num_extbc = arglist.num_extbc;
		ExtBC::extBC<T> extbc[num_extbc];
		vector<int> *relative_bc_index = new vector<int> [size];
		
		if (arglist.num_extbc > 0)
		{
			for(int i=0; i<arglist.num_extbc; i++){
				extbc[i] = ExtBC::extBC<T>(arglist.extbc_fname[i], arglist.extbc_bctype[i]);
				if(arglist.extbc_bctype[i]==1){
					extbc[i].convert_to_secs();
				}
			}
			
			for(int i=0; i<arglist.num_extbc; i++){
				extbc[i].extreme_cols.assign(2, 0); //there are two points per extbc
				extbc[i].extreme_rows.assign(2, 0); //there are two points per extbc
				
				T extreme_x1 = arglist.extbc_x1_loc[i];
				T extreme_y1 = arglist.extbc_y1_loc[i];
				T extreme_x2 = arglist.extbc_x2_loc[i];
				T extreme_y2 = arglist.extbc_y2_loc[i];

				extbc[i].extreme_cols[0] = calc_src_col(extreme_x1, dem.get_xll_corner(), dem.get_cell_size());
				extbc[i].extreme_rows[0] = calc_src_row(extreme_y1, dem.get_yll_corner(), dem.get_cell_size(), org_rows);
				extbc[i].extreme_cols[1] = calc_src_col(extreme_x2, dem.get_xll_corner(), dem.get_cell_size());
				extbc[i].extreme_rows[1] = calc_src_row(extreme_y2, dem.get_yll_corner(), dem.get_cell_size(), org_rows);	

				extbc[i].ncells=extbc[i].check_extreme_extbc(extbc[i].extreme_cols,extbc[i].extreme_rows,org_cols,org_rows);
				extbc[i].ncells_local=0;
			}
			
			for(int i=0; i<arglist.num_extbc; i++){
				extbc[i].create_involved_cells(extbc[i].extreme_cols,extbc[i].extreme_rows,org_cols,org_rows,arglist.extbc_bctype[i]);
				if(strcmp(arglist.input_option.c_str(), "SEQ")==0){
					if(rank==0){
						dem.copy_elevation_into_ghost_cells(extbc[i].i_rows,extbc[i].i_cols,extbc[i].ncells, extbc[i].location);
					}
				}else{
					//in parallel input mode, we cannot copy the elevation into ghost cells because the sub_dem files haven't been loaded yet. It's done after while reading in parallel				
				}
			}
			
			for(int i=0; i<arglist.num_extbc; i++){
				for(int j=0; j<extbc[i].ncells; j++){
					int srank = 0;
					int prev_rows_sum = 0;
					int cell_row = extbc[i].i_rows[j];
					if(size > 1){
						int rows_sum = pd.part_dims[0].first - 2 * GHOST_CELL_PADDING;
						if(cell_row >= rows_sum){
							for(int k=1; k<size; k++){
								prev_rows_sum = rows_sum;
								rows_sum += pd.part_dims[k].first - 2 * GHOST_CELL_PADDING;
								if(cell_row < rows_sum){
									srank = k;
									break;
								}
							}
						}
					}
					cell_row = cell_row - prev_rows_sum + GHOST_CELL_PADDING;
					int new_index = cell_row * (org_cols+2*GHOST_CELL_PADDING) + extbc[i].i_cols[j]+GHOST_CELL_PADDING;
					relative_bc_index[srank].push_back(new_index);
					if (rank == srank){
						extbc[i].ncells_local++;
					}

				}
			}
		}
		
		num_of_extbc = arglist.num_extbc;
		num_extbc_cells = relative_bc_index[rank].size();
		
		host_bc_cells_size = 0;
		if (num_of_extbc > 0)
		{
			for(int i=0; i<num_of_extbc; i++){
				host_bc_cells_size += extbc[i].ncells_local;
			}
		}else{
			host_bc_cells_size = 1;
		}

		host_relative_bc_index = new int[max(host_bc_cells_size,1)];
		host_bc_type = new int[max(host_bc_cells_size,1)];
		host_bc_start_index = new int[max(host_bc_cells_size,1)];
		host_bc_nrows_vars = new int[max(host_bc_cells_size,1)];

		if (num_of_extbc > 0)
		{
			int moving_index = 0;
			int start_index = 0;
			for(int i=0; i<num_of_extbc; i++){
				int extbc_num_rows=extbc[i].get_num_rows();
				for(int j=0; j<extbc[i].ncells_local; j++){
					host_relative_bc_index[moving_index+j] = relative_bc_index[rank][moving_index+j];
					host_bc_type[moving_index+j] = arglist.extbc_bctype[i];
					host_bc_start_index[moving_index+j] = start_index;
					host_bc_nrows_vars[moving_index+j] = extbc_num_rows;
				}
				start_index +=extbc_num_rows;
				moving_index += extbc[i].ncells_local;
			}
		}
		
		host_bc_vars_arr_size = 0;
		if (num_of_extbc > 0)
		{
			for(int i=0; i<num_of_extbc; i++){
				host_bc_vars_arr_size += extbc[i].get_num_rows();
			}
		}else{
			host_bc_vars_arr_size = 1;
		}


		host_extbc_var1_arr = new T[host_bc_vars_arr_size];
		host_extbc_var2_arr = new T[host_bc_vars_arr_size];

		if (num_of_extbc > 0)
		{
			int moving_index2 = 0;
			for(int i=0; i<num_of_extbc; i++){
				int extbc_num_rows=extbc[i].get_num_rows();
				for (int j = 0; j < extbc_num_rows; j++)
				{
					host_extbc_var1_arr[moving_index2+j] = extbc[i].get_var1_at(j);
					host_extbc_var2_arr[moving_index2+j] = extbc[i].get_var2_at(j);
				}
				moving_index2 += extbc_num_rows;
			}

		}else{
			host_extbc_var1_arr[0]=0.0;
			host_extbc_var2_arr[0]=0.0;
		}
		
		delete[] relative_bc_index;
	}


	template<typename T>
	void triton<T>::partition_matrix_files()
	{

		//Note that the first partitioning at the beginning of the simulation is homogeneous (static)	
		if(size > 1 && rank == 0){
			std::cerr << IN "Creating partition data" << std::endl;
		}

		if(size > 1)
		{
			sub_dem = MpiUtils::scatter_exchange(dem.get_data(), pd, rank);
			sub_nin = MpiUtils::scatter_exchange(nin.get_data(), pd, rank);
		}
		else
		{
			sub_dem = dem;
			sub_nin = nin;
		}
		
		rows = sub_dem.get_num_rows();
		cols = sub_dem.get_num_cols();
		
		sub_dem.set_nrows(sub_dem.get_num_rows());
		sub_dem.set_ncols(sub_dem.get_num_cols());
		sub_dem.set_cell_size(dem.get_cell_size());
		sub_dem.set_xll_corner(dem.get_xll_corner());
		sub_dem.set_yll_corner(dem.get_yll_corner());
		sub_dem.set_no_data_value(dem.get_no_data_value());
		
		if(arglist.runoff_map.size() > 0)
		{
			if(size > 1)
			{
				sub_rin = MpiUtils::scatter_exchange_int(rin.get_data(), pd, rank);
			}
			else
			{
				sub_rin = rin;
			}
		}
		
		if(arglist.h_infile.size() > 0)
		{
			if(size > 1)
			{
				sub_hin = MpiUtils::scatter_exchange(hin.get_data(), pd, rank);
			}
			else
			{
				sub_hin = hin;
			}
		}
		
		if(arglist.qx_infile.size() > 0)
		{
			if(size > 1)
			{
				sub_qxin = MpiUtils::scatter_exchange(uin.get_data(), pd, rank);
			}
			else
			{
				sub_qxin = uin;
			}
		}
		
		if(arglist.qy_infile.size() > 0)
		{
			if(size > 1)
			{
				sub_qyin = MpiUtils::scatter_exchange(vin.get_data(), pd, rank);
			}
			else
			{
				sub_qyin = vin;
			}
		}
	
		if (arglist.checkpoint_id > 0)
		{
			if(strcmp(arglist.output_option.c_str(), "SEQ") == 0)
			{
				if(size > 1)
				{
					sub_hot_hin = MpiUtils::scatter_exchange(hot_hin.get_data(), pd, rank);
					sub_hot_qxin = MpiUtils::scatter_exchange(hot_qxin.get_data(), pd, rank);
					sub_hot_qyin = MpiUtils::scatter_exchange(hot_qyin.get_data(), pd, rank);
					
					if (arglist.max_value_print_option.size() > 0)	
					{
						sub_max_value_h = MpiUtils::scatter_exchange(hot_max_value_h.get_data(), pd, rank);
					}
				}
				else
				{
					sub_hot_hin = hot_hin;
					sub_hot_qxin = hot_qxin;
					sub_hot_qyin = hot_qyin;
					
					if (arglist.max_value_print_option.size() > 0)	
					{
						sub_max_value_h = hot_max_value_h;
					}
				}
			}
			else
			{
				string temp_num(to_string(arglist.checkpoint_id));
				if (arglist.checkpoint_id < 10)
				{
					temp_num = "0" + temp_num;
				}
				
				string temp_rank(to_string(rank));
				if (rank < 10)
				{
					temp_rank = "0" + temp_rank;
				}

				int host_dem_original_row = rows - 2 * GHOST_CELL_PADDING;
				int host_dem_original_col = cols - 2 * GHOST_CELL_PADDING;

				string filedirH(project_dir + "/" + OUTPUT_DIR + "/" + BIN_DIR + "/H_" + temp_num + "_" + temp_rank + ".out");
				sub_hot_hin.load_from_binary_file(host_dem_original_row, host_dem_original_col, filedirH);
				sub_hot_hin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
				if(arglist.open_boundaries){				
					sub_hot_hin.copy_value_into_ghost_cells();
				}
				string filedirU(project_dir + "/" + OUTPUT_DIR + "/" + BIN_DIR + "/QX_" + temp_num + "_" + temp_rank + ".out");
				sub_hot_qxin.load_from_binary_file(host_dem_original_row, host_dem_original_col, filedirU);
				sub_hot_qxin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
				if(arglist.open_boundaries){				
					sub_hot_qxin.copy_value_into_ghost_cells();
				}
				string filedirV(project_dir + "/" + OUTPUT_DIR + "/" + BIN_DIR + "/QY_" + temp_num + "_" + temp_rank + ".out");
				sub_hot_qyin.load_from_binary_file(host_dem_original_row, host_dem_original_col, filedirV);
				sub_hot_qyin.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
				if(arglist.open_boundaries){				
					sub_hot_qyin.copy_value_into_ghost_cells();
				}
				MpiUtils::exchange(sub_hot_hin.begin(), rows, cols, rank, size, USE_MATRIX);
				MPI_Barrier(MPI_COMM_WORLD);
				MpiUtils::exchange(sub_hot_qxin.begin(), rows, cols, rank, size, USE_MATRIX);
				MPI_Barrier(MPI_COMM_WORLD);
				MpiUtils::exchange(sub_hot_qyin.begin(), rows, cols, rank, size, USE_MATRIX);
				MPI_Barrier(MPI_COMM_WORLD);
				
				if (arglist.max_value_print_option.size() > 0)	
				{
					string filedirMaxH(project_dir + "/" + OUTPUT_DIR + "/" + BIN_DIR + "/MH_" + temp_num + "_" + temp_rank + ".out");
					sub_max_value_h.load_from_binary_file(host_dem_original_row, host_dem_original_col, filedirMaxH);
					sub_max_value_h.add_ghost_cells(GHOST_CELL_PADDING, GHOST_CELL_PADDING, 0.0);
					if(arglist.open_boundaries){				
						sub_max_value_h.copy_value_into_ghost_cells();
					}
				}
			}
		}

		//overwrite hin, qxin, qyin 
		if (arglist.checkpoint_id > 0)
		{
			sub_hin = sub_hot_hin;
			sub_qxin = sub_hot_qxin;
			sub_qyin = sub_hot_qyin;
			
			if (arglist.max_value_print_option.size() <= 0)	
			{
				sub_max_value_h = sub_hin;
			}
		}
		else
		{
			Matrix::matrix<T> values(rows, cols);
			values.zero_fill();

			if (arglist.h_infile.size() <= 0)
			{
				sub_hin = values;
			}
			if (arglist.qx_infile.size() <= 0)
			{
				sub_qxin = values;
			}
			if (arglist.qy_infile.size() <= 0)
			{
				sub_qyin = values;
			}
			
			sub_max_value_h = sub_hin;
		}
		if(size > 1 && rank == 0){
			std::cerr << OK "Data has been partitioned: " << arglist.domain_decomposition << " domain decomposition" <<std::endl;
		}

	}


	template<typename T>
	void triton<T>::process_runoff()
	{
		if (arglist.num_runoffs > 0)
		{
			index_row_runoff = 0;
			for (int j = 0; j < roff.get_num_inflow_rows() - 1; j++)
			{
				if (simtime > roff.get_time_at(j) && simtime <= roff.get_time_at(j + 1))
				{
					index_row_runoff = j;
				}
			}
		}

		int num_runoffs_temp = arglist.num_runoffs;
		if (num_runoffs_temp == 0)
		{
			num_runoffs_temp = 1;
		}
		
		int runoff_row_size_temp = 1;
		if (arglist.num_runoffs > 0)
		{
			runoff_row_size_temp = roff.get_num_inflow_rows();
		}
		host_runoff_id_arr = new int[rows * cols];
		host_runoff_intensity_arr_size = num_runoffs_temp * runoff_row_size_temp;
		host_runoff_intensity_arr = new T[host_runoff_intensity_arr_size];

		for (int j = 0; j < rows * cols; j++)
		{
			if (arglist.num_runoffs > 0)
			{
				host_runoff_id_arr[j] = sub_rin.get_value(j);
			}
			else
			{
				host_runoff_id_arr[j] = -1;
			}
		}

		for (int j = 0; j < num_runoffs_temp * runoff_row_size_temp; j++)
		{
			host_runoff_intensity_arr[j] = 0.0;
		}
		
		
		if (arglist.num_runoffs > 0)
		{
			for (int j = 0; j < arglist.num_runoffs; j++)
			{
				for (int k = 0; k < roff.get_num_inflow_rows(); k++)
				{
					host_runoff_intensity_arr[j*roff.get_num_inflow_rows() + k] = roff.get_flow_at(k, j + 1);
				}
			}
		}					

	}


	template<typename T>
	void triton<T>::create_host_aux_vectors()
	{
		host_halo_arr_size = 12 * cols*GHOST_CELL_PADDING;
		host_halo_arr = new T[host_halo_arr_size];
		for (int i = 0; i < host_halo_arr_size; i++)
		{
			host_halo_arr[i] = 0.0;
		}

		host_sqrth_arr = new T[rows * cols];
		for (int j = 0; j < rows * cols; j++)
		{
			host_sqrth_arr[j] = 0.0;
		}
		
#ifdef ACTIVE_GPU
		if ((rows*cols) % THREAD_BLOCK == 0)
		{
			host_reduce_dt_arr_sz = (rows*cols) / THREAD_BLOCK;
		}
		else
		{
			host_reduce_dt_arr_sz = (rows*cols) / THREAD_BLOCK + 1;
		}
#else
		host_reduce_dt_arr_sz = rows * cols;
#endif
		
		host_dt_values_arr = new T[host_reduce_dt_arr_sz];
		for (int i = 0; i < host_reduce_dt_arr_sz; i++)
		{
			host_dt_values_arr[i] = 0.0;
		}
		
		
		host_rhsh0 = new T[rows*cols];
		host_rhsh1 = new T[rows*cols];
		host_rhsqx0 = new T[rows*cols];
		host_rhsqx1 = new T[rows*cols];
		host_rhsqy0 = new T[rows*cols];
		host_rhsqy1 = new T[rows*cols];

		for (int i = 0; i < rows*cols; i++)
		{
			host_rhsh0[i] = 0.0;
			host_rhsh1[i] = 0.0;
			host_rhsqx0[i] = 0.0;
			host_rhsqx1[i] = 0.0;
			host_rhsqy0[i] = 0.0;
			host_rhsqy1[i] = 0.0;
		}
	}


	template<typename T>
	void triton<T>::create_host_vectors()
	{
		host_vec = std::vector<T*>();
		host_vec_int = std::vector<int*>();

		host_vec.push_back(sub_hin.get_data());
		host_vec.push_back(sub_qxin.get_data());
		host_vec.push_back(sub_qyin.get_data());
		host_vec.push_back(sub_nin.get_data());
		host_vec.push_back(sub_dem.get_data());
		host_vec.push_back(sub_max_value_h.get_data());

		host_vec.push_back(host_rhsh0);
		host_vec.push_back(host_rhsh1);
		host_vec.push_back(host_rhsqx0);
		host_vec.push_back(host_rhsqx1);
		host_vec.push_back(host_rhsqy0);
		host_vec.push_back(host_rhsqy1);

		host_vec.push_back(host_sqrth_arr);
		host_vec.push_back(host_halo_arr);
		host_vec.push_back(host_dt_values_arr);
		host_vec.push_back(host_hyg_time_arr);
		host_vec.push_back(host_hyg_val_arr);
		host_vec.push_back(host_runoff_intensity_arr);
		host_vec.push_back(host_extbc_var1_arr);
		host_vec.push_back(host_extbc_var2_arr);

		host_vec_int.push_back(host_src_pos_arr);
		host_vec_int.push_back(host_runoff_id_arr);
		host_vec_int.push_back(host_relative_bc_index);
		host_vec_int.push_back(host_bc_type);
		host_vec_int.push_back(host_bc_start_index);
		host_vec_int.push_back(host_bc_nrows_vars);
	}


	template<typename T>
	void triton<T>::create_device_vectors()
	{
		nbytes = (sizeof(T) * rows * cols);
		nbytes_halo = (sizeof(T) * host_halo_arr_size);

#ifdef ACTIVE_GPU

		int deviceId = 0;
		cudaError_t err = cudaGetDevice(&deviceId);
		if (err != cudaSuccess) 
		{
			std::cerr << cudaGetErrorString(err) << std::endl;
			exit(EXIT_FAILURE);
		}
		
		int deviceCount = 0;
		err = cudaGetDeviceCount(&deviceCount);
		if (err != cudaSuccess) 
		{
			std::cerr << cudaGetErrorString(err) << std::endl;
			exit(EXIT_FAILURE);
		}
		
		if(deviceId != rank % deviceCount)
		{
			err = cudaSetDevice(rank % deviceCount);
			if (err != cudaSuccess) 
			{
				std::cerr << cudaGetErrorString(err) << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		cudaStreamCreate(&streams);

		device_vec = std::vector<T*>();
		device_vec_int = std::vector<int*>();

		int nbytes_dt = (sizeof(T) * host_reduce_dt_arr_sz);
		int nbytes_hyg_time = (sizeof(T) * host_hyg_time_arr_size);
		int nbytes_hyg_val = (sizeof(T) * host_hyg_time_arr_size * host_src_pos_arr_size);
		int nbytes_runoff_intensity = (sizeof(T) * host_runoff_intensity_arr_size);
		int nbytes_src_pos = (sizeof(int) * host_src_pos_arr_size);
		int nbytes_runoff_id = (sizeof(int) * rows * cols);
		int nbytes_bc_cell_size = (sizeof(int) * max(host_bc_cells_size,1));
		int nbytes_bc_vars = (sizeof(T) * host_bc_vars_arr_size);


		T *device_h, *device_qx, *device_qy,
		*device_n, *device_dem, *device_sqrth_arr, *device_halo_arr, *device_dt_values_arr,
		*device_rhsh0, *device_rhsh1, *device_rhsqx0, *device_rhsqx1, *device_rhsqy0, *device_rhsqy1,
		*device_hyg_time_arr, *device_hyg_val_arr,
		*device_runoff_intensity_arr,
		*device_bc_var1_arr, *device_bc_var2_arr, *device_max_value_h;

		int *device_src_pos_arr, *device_runoff_id_arr, *device_relative_bc_index, *device_bc_type,
		*device_bc_start_index, *device_bc_nrows_vars;
		
		cudaMalloc((void**)&device_h, nbytes);
		cudaMalloc((void**)&device_qx, nbytes);
		cudaMalloc((void**)&device_qy, nbytes);
		cudaMalloc((void**)&device_n, nbytes);
		cudaMalloc((void**)&device_dem, nbytes);
		cudaMalloc((void**)&device_max_value_h, nbytes);

		cudaMalloc((void**)&device_rhsh0, nbytes);
		cudaMalloc((void**)&device_rhsh1, nbytes);
		cudaMalloc((void**)&device_rhsqx0, nbytes);
		cudaMalloc((void**)&device_rhsqx1, nbytes);
		cudaMalloc((void**)&device_rhsqy0, nbytes);
		cudaMalloc((void**)&device_rhsqy1, nbytes);

		cudaMalloc((void**)&device_sqrth_arr, nbytes);
		cudaMalloc((void**)&device_halo_arr, nbytes_halo);
		cudaMalloc((void**)&device_dt_values_arr, nbytes_dt);
		cudaMalloc((void**)&device_hyg_time_arr, nbytes_hyg_time);
		cudaMalloc((void**)&device_hyg_val_arr, nbytes_hyg_val);
		cudaMalloc((void**)&device_runoff_intensity_arr, nbytes_runoff_intensity);
		cudaMalloc((void**)&device_bc_var1_arr, nbytes_bc_vars);
		cudaMalloc((void**)&device_bc_var2_arr, nbytes_bc_vars);


		cudaMalloc((void**)&device_src_pos_arr, nbytes_src_pos);
		cudaMalloc((void**)&device_runoff_id_arr, nbytes_runoff_id);
		cudaMalloc((void**)&device_relative_bc_index, nbytes_bc_cell_size);
		cudaMalloc((void**)&device_bc_type, nbytes_bc_cell_size);
		cudaMalloc((void**)&device_bc_start_index, nbytes_bc_cell_size);
		cudaMalloc((void**)&device_bc_nrows_vars, nbytes_bc_cell_size);

		device_vec.push_back(device_h);
		device_vec.push_back(device_qx);
		device_vec.push_back(device_qy);
		device_vec.push_back(device_n);
		device_vec.push_back(device_dem);
		device_vec.push_back(device_max_value_h);

		device_vec.push_back(device_rhsh0);
		device_vec.push_back(device_rhsh1);
		device_vec.push_back(device_rhsqx0);
		device_vec.push_back(device_rhsqx1);
		device_vec.push_back(device_rhsqy0);
		device_vec.push_back(device_rhsqy1);

		device_vec.push_back(device_sqrth_arr);
		device_vec.push_back(device_halo_arr);
		device_vec.push_back(device_dt_values_arr);
		device_vec.push_back(device_hyg_time_arr);
		device_vec.push_back(device_hyg_val_arr);
		device_vec.push_back(device_runoff_intensity_arr);
		device_vec.push_back(device_bc_var1_arr);
		device_vec.push_back(device_bc_var2_arr);

		device_vec_int.push_back(device_src_pos_arr);
		device_vec_int.push_back(device_runoff_id_arr);
		device_vec_int.push_back(device_relative_bc_index);
		device_vec_int.push_back(device_bc_type);
		device_vec_int.push_back(device_bc_start_index);
		device_vec_int.push_back(device_bc_nrows_vars);


		cudaMemcpyAsync(device_vec[H], host_vec[H], nbytes, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[QX], host_vec[QX], nbytes, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[QY], host_vec[QY], nbytes, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[NMAN], host_vec[NMAN], nbytes, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[DEM], host_vec[DEM], nbytes, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[MAXH], host_vec[MAXH], nbytes, cudaMemcpyHostToDevice, streams);

		cudaMemcpyAsync(device_vec[RHSH0], host_vec[RHSH0], nbytes, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[RHSH1], host_vec[RHSH1], nbytes, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[RHSQX0], host_vec[RHSQX0], nbytes, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[RHSQX1], host_vec[RHSQX1], nbytes, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[RHSQY0], host_vec[RHSQY0], nbytes, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[RHSQY1], host_vec[RHSQY1], nbytes, cudaMemcpyHostToDevice, streams);

		cudaMemcpyAsync(device_vec[SQRTH], host_vec[SQRTH], nbytes, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[HALO], host_vec[HALO], nbytes_halo, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[DT], host_vec[DT], nbytes_dt, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[HYGT], host_vec[HYGT], nbytes_hyg_time, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[HYGV], host_vec[HYGV], nbytes_hyg_val, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[RUNIN], host_vec[RUNIN], nbytes_runoff_intensity, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[EXTBCV1], host_vec[EXTBCV1], nbytes_bc_vars, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec[EXTBCV2], host_vec[EXTBCV2], nbytes_bc_vars, cudaMemcpyHostToDevice, streams);


		cudaMemcpyAsync(device_vec_int[SRCP], host_vec_int[SRCP], nbytes_src_pos, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec_int[RUNID], host_vec_int[RUNID], nbytes_runoff_id, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec_int[BCRELATIVEINDEX], host_vec_int[BCRELATIVEINDEX], nbytes_bc_cell_size, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec_int[BCTYPE], host_vec_int[BCTYPE], nbytes_bc_cell_size, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec_int[BCINDEXSTART], host_vec_int[BCINDEXSTART], nbytes_bc_cell_size, cudaMemcpyHostToDevice, streams);
		cudaMemcpyAsync(device_vec_int[BCNROWSVARS], host_vec_int[BCNROWSVARS], nbytes_bc_cell_size, cudaMemcpyHostToDevice, streams);
		cudaStreamSynchronize(streams);
#endif
	}


	template<typename T>
	int triton<T>::calc_src_col(T src_x, T xllc, T cell_size_)
	{
		return ceil(((src_x - xllc) / cell_size_)+1e-16) - 1;
	}


	template<typename T>
	int triton<T>::calc_src_row(T src_y, T yllc, T cell_size_, int nrows)
	{
		return ceil((nrows - ((src_y - yllc) / cell_size_))+1e-16) - 1;
	}


	template<class T>
	triton<T>::~triton()
	{		
		//dont delete matrix object
		delete[] host_vec[EXTBCV2];
		delete[] host_vec[EXTBCV1];
		delete[] host_vec[RUNIN];
		delete[] host_vec[HYGV];
		delete[] host_vec[HYGT];
		delete[] host_vec[DT];
		delete[] host_vec[HALO];
		delete[] host_vec[SQRTH];
		delete[] host_vec[RHSQY1];
		delete[] host_vec[RHSQY0];
		delete[] host_vec[RHSQX1];
		delete[] host_vec[RHSQX0];
		delete[] host_vec[RHSH1];
		delete[] host_vec[RHSH0];

		delete[] host_vec_int[BCNROWSVARS];
		delete[] host_vec_int[BCINDEXSTART];
		delete[] host_vec_int[BCTYPE];
		delete[] host_vec_int[BCRELATIVEINDEX];
		delete[] host_vec_int[RUNID];
		delete[] host_vec_int[SRCP];

#ifdef ACTIVE_GPU
		while (!device_vec.empty())
		{
			cudaFree(device_vec.back());
			device_vec.pop_back();
		}

		while (!device_vec_int.empty())
		{
			cudaFree(device_vec_int.back());
			device_vec_int.pop_back();
		}
#endif
	}


	template<typename T>
	void triton<T>::simulate()
	{
		compute_init_dt();
		
		if(rank==0){
			std::cerr << OK "Simulation starts" << std::endl;
		}

		MPI_Barrier(MPI_COMM_WORLD);


		st.start(SIMULATION_TIME);

		out.init(rows, cols, rank, size, project_dir, arglist.outfile_pattern, arglist.time_series_flag, cfg_content, arglist.output_option);

		if (arglist.time_series_flag)
		{
			out.init_time_series(num_of_obs_points, arglist.observation_x_loc.size(), relative_obs_index, observation_cells, observation_cells_global);
		}


		global_dt = arglist.time_step;
		average_dt = 0.0;
		local_dt = arglist.time_step;
		int it_count = arglist.it_count;
		int it_count_average = 0;
		int print_id = arglist.checkpoint_id;
		
		//this is to allow simtime different from zero without checkpointing
		if(print_id==0 && simtime>0.0){
			print_id=simtime/arglist.print_interval;
		}
		
		while (simtime < arglist.sim_duration)
		{
			it_count++;
			it_count_average++;

			if (!arglist.time_increment_fixed)
			{
				compute_local_dt();
				compute_global_dt(print_id);
			}
			compute_new_state();

			simtime += global_dt;
			average_dt+=global_dt;

			if (simtime >= arglist.print_interval * (print_id + 1))
			{
				print_id++;

#ifdef ACTIVE_GPU
				st.start(COMPUTE_TIME);
				cudaMemcpyAsync(host_vec[H], device_vec[H], nbytes, cudaMemcpyDeviceToHost, streams);
				cudaMemcpyAsync(host_vec[QX], device_vec[QX], nbytes, cudaMemcpyDeviceToHost, streams);
				cudaMemcpyAsync(host_vec[QY], device_vec[QY], nbytes, cudaMemcpyDeviceToHost, streams);
				if (arglist.max_value_print_option.size() > 0)
				{
					cudaMemcpyAsync(host_vec[MAXH], device_vec[MAXH], nbytes, cudaMemcpyDeviceToHost, streams);
				}
				cudaStreamSynchronize(streams);
				st.stop(COMPUTE_TIME);
#endif

				st.start(IO_TIME);
				out.write_output(sub_hin, sub_qxin, sub_qyin, arglist.output_format, arglist.print_option, print_id, it_count, simtime, average_dt/it_count_average,sub_max_value_h, arglist.max_value_print_option);
				it_count_average=0;
				average_dt=0.0;
				st.stop(IO_TIME);

				#if WRITE_PERFORMANCE
					st.stop(SIMULATION_TIME);
					st.stop(TOTAL_TIME);
					out.write_times(st, print_id);
					st.start(SIMULATION_TIME);
					st.start(TOTAL_TIME);
				#endif
				
				//mandatory in case of hotstart
				if(strcmp(arglist.domain_decomposition.c_str(), TYPE_DYNAMIC)==0 && size > 1){
					out.write_domain_decomposition(pd,print_id);
				}

				//we want to avoid repartitioning in the last iteration when simtime=arglist.sim_duration
				if(strcmp(arglist.domain_decomposition.c_str(), TYPE_DYNAMIC)==0 && print_id%arglist.factor_interval_domain_decomposition==0 && size > 1 && simtime < arglist.sim_duration)
				{
					st.start(RESIZE_TIME);
					new_domain_decomposition();
					st.stop(RESIZE_TIME);
				}

			}

		}
		st.stop(SIMULATION_TIME);
		st.stop(TOTAL_TIME);
		
		out.write_times(st, -1);
		if(rank==0){
			std::cerr << OK "Simulation ends" << std::endl;
		}


	}


	template<typename T>
	void triton<T>::compute_init_dt()
	{
		init_dt=arglist.print_interval;
		if (arglist.num_runoffs > 0){
			init_dt=fmin(init_dt,roff.get_time_at(1)-roff.get_time_at(0));
		}
		if(num_of_src>0){
			init_dt=fmin(init_dt,hyg.get_time_at(1)-hyg.get_time_at(0));
		}
		init_dt*=0.01;

	}


	template<typename T>
	void triton<T>::compute_local_dt()
	{
		st.start(COMPUTE_TIME);
#ifdef ACTIVE_GPU
		int cur_dt_arr_sz = host_reduce_dt_arr_sz;

		Kernels::compute_dt_and_sqrt << <(rows*cols + THREAD_BLOCK - 1) / THREAD_BLOCK, THREAD_BLOCK, THREAD_BLOCK * sizeof(T), streams >> > (rows*cols, cell_size,
		device_vec[QX], device_vec[QY], device_vec[H], device_vec[SQRTH], device_vec[DT], arglist.courant, arglist.hextra);

		while (cur_dt_arr_sz > 1)
		{
			int temp_dt_arr_sz = (cur_dt_arr_sz / THREAD_BLOCK) + 1;
			if(cur_dt_arr_sz % THREAD_BLOCK == 0)
			{
				temp_dt_arr_sz = (cur_dt_arr_sz / THREAD_BLOCK);
			}
			Kernels::find_min_dt << <temp_dt_arr_sz, THREAD_BLOCK, THREAD_BLOCK * sizeof(T), streams >> > (cur_dt_arr_sz, device_vec[DT]);

			cur_dt_arr_sz = temp_dt_arr_sz;
		}

		cudaMemcpyAsync(&local_dt, device_vec[DT], sizeof(T), cudaMemcpyDeviceToHost, streams);
		cudaStreamSynchronize(streams);
#else
		Kernels::compute_dt_and_sqrt(rows*cols, cell_size, host_vec[QX], host_vec[QY], host_vec[H], host_vec[SQRTH], host_vec[DT], arglist.courant, arglist.hextra);
		Kernels::find_min_dt(rows*cols, host_vec[DT]);
		local_dt = host_vec[DT][0];
#endif
		st.stop(COMPUTE_TIME);
	}


	template<typename T>
	void triton<T>::compute_global_dt(int print_id)
	{
		if (size > 1)
		{
			st.start(BALANCING_MPI_TIME);
			st.start(MPI_TIME);
			MPI_Allreduce(&local_dt, &global_dt, 1, MPI_DATA_TYPE, MPI_MIN, MPI_COMM_WORLD);
			st.stop(MPI_TIME);
			st.stop(BALANCING_MPI_TIME);
		}
		else
		{
			global_dt = local_dt;
		}

		if (global_dt >= MAX_VALUE - 1.0)
		{
			global_dt = init_dt;
		}

		if (simtime + global_dt > arglist.print_interval * (print_id + 1))
		{
			global_dt = arglist.print_interval * (print_id + 1) - simtime;
		}
	}


	template<typename T>
	void triton<T>::compute_new_state()
	{
		st.start(COMPUTE_TIME);

#ifdef ACTIVE_GPU
		Kernels::flux_x << <(rows*cols + THREAD_BLOCK - 1) / THREAD_BLOCK, THREAD_BLOCK, 0, streams >> > (rows*cols, rows, cols, cell_size, global_dt,
		device_vec[H], device_vec[QX], device_vec[QY], device_vec[DEM], device_vec[SQRTH],
		device_vec[RHSH0], device_vec[RHSH1], device_vec[RHSQX0], device_vec[RHSQX1], device_vec[RHSQY0], device_vec[RHSQY1], arglist.hextra);

		Kernels::flux_y << <(rows*cols + THREAD_BLOCK - 1) / THREAD_BLOCK, THREAD_BLOCK, 0, streams >> > (rows*cols, rows, cols, cell_size, global_dt,
		device_vec[H], device_vec[QX], device_vec[QY], device_vec[DEM], device_vec[SQRTH],
		device_vec[RHSH0], device_vec[RHSH1], device_vec[RHSQX0], device_vec[RHSQX1], device_vec[RHSQY0], device_vec[RHSQY1], arglist.hextra);

		Kernels::update_cells << <(rows*cols + THREAD_BLOCK - 1) / THREAD_BLOCK, THREAD_BLOCK, 0, streams >> > (rows*cols, rows, cols, global_dt,
		device_vec[H], device_vec[QX], device_vec[QY], device_vec[DEM], device_vec[NMAN],
		device_vec[RHSH0], device_vec[RHSH1], device_vec[RHSQX0], device_vec[RHSQX1], device_vec[RHSQY0], device_vec[RHSQY1], arglist.hextra);
#else
		Kernels::flux_x(rows*cols, rows, cols, cell_size, global_dt,
		host_vec[H], host_vec[QX], host_vec[QY], host_vec[DEM], host_vec[SQRTH],
		host_vec[RHSH0], host_vec[RHSH1], host_vec[RHSQX0], host_vec[RHSQX1], host_vec[RHSQY0], host_vec[RHSQY1], arglist.hextra);

		Kernels::flux_y(rows*cols, rows, cols, cell_size, global_dt,
		host_vec[H], host_vec[QX], host_vec[QY], host_vec[DEM], host_vec[SQRTH],
		host_vec[RHSH0], host_vec[RHSH1], host_vec[RHSQX0], host_vec[RHSQX1], host_vec[RHSQY0], host_vec[RHSQY1], arglist.hextra);

		Kernels::update_cells(rows*cols, rows, cols, global_dt,
		host_vec[H], host_vec[QX], host_vec[QY], host_vec[DEM], host_vec[NMAN],
		host_vec[RHSH0], host_vec[RHSH1], host_vec[RHSQX0], host_vec[RHSQX1], host_vec[RHSQY0], host_vec[RHSQY1], arglist.hextra);
#endif

		if (arglist.num_runoffs > 0)
		{
			if(simtime > roff.get_time_at(roff.get_num_inflow_rows()-1))
			{	
				index_row_runoff=roff.get_num_inflow_rows()-1;
			}
			else
			{
				if (simtime > roff.get_time_at(index_row_runoff + 1))
				{
					index_row_runoff++;
				}
			}

#ifdef ACTIVE_GPU
			Kernels::update_runoff << <(rows*cols + THREAD_BLOCK - 1) / THREAD_BLOCK, THREAD_BLOCK, 0, streams >> > (rows*cols, rows, cols, global_dt,
			device_vec_int[RUNID], index_row_runoff, roff.get_num_inflow_rows(), device_vec[RUNIN], device_vec[H], device_vec[QX], device_vec[QY], arglist.hextra);
#else
			Kernels::update_runoff(rows*cols, rows, cols, global_dt,
			host_vec_int[RUNID], index_row_runoff, roff.get_num_inflow_rows(), host_vec[RUNIN], host_vec[H], host_vec[QX], host_vec[QY], arglist.hextra);
#endif
		}

		if (num_of_src > 0)
		{
			bool check_flag = false;
			int idx_high = hyg.get_num_inflow_rows() - 1;
			for (int i = idx_low; i < idx_high; i++)
			{
				if (hyg.get_time_at(i + 1) > simtime)
				{
					idx_low = i;
					idx_high = i + 1;
					check_flag = true;
					break;
				}
			}
			if (!check_flag)
			{
				idx_low = idx_high;
			}

#ifdef ACTIVE_GPU
			Kernels::compute_flow << <(num_of_src + THREAD_BLOCK - 1) / THREAD_BLOCK, THREAD_BLOCK, 0, streams >> > (num_of_src, device_vec[HYGT], device_vec[HYGV],
			cell_size, global_dt, simtime, idx_low, idx_high, device_vec[H], device_vec_int[SRCP]);
#else
			Kernels::compute_flow(num_of_src, host_vec[HYGT], host_vec[HYGV],
			cell_size, global_dt, simtime, idx_low, idx_high, host_vec[H], host_vec_int[SRCP]);
#endif
		}


		if(arglist.open_boundaries){

			#ifdef ACTIVE_GPU
				Kernels::copy_info_to_exterior_boundaries_west_east << <(2*rows*GHOST_CELL_PADDING + THREAD_BLOCK - 1) / THREAD_BLOCK, THREAD_BLOCK, 0, streams >> > (2*rows*GHOST_CELL_PADDING, rows, cols, device_vec[H], device_vec[QX], device_vec[QY]);
				Kernels::copy_info_to_exterior_boundaries_north_south << <(cols*GHOST_CELL_PADDING + THREAD_BLOCK - 1) / THREAD_BLOCK, THREAD_BLOCK, 0, streams >> > (cols*GHOST_CELL_PADDING, rows, cols, device_vec[H], device_vec[QX], device_vec[QY],rank, size);
			#else
				Kernels::copy_info_to_exterior_boundaries_west_east(2*rows*GHOST_CELL_PADDING, rows, cols, host_vec[H], host_vec[QX], host_vec[QY]);
				Kernels::copy_info_to_exterior_boundaries_north_south(cols*GHOST_CELL_PADDING, rows, cols, host_vec[H], host_vec[QX], host_vec[QY],rank, size);
			#endif

		}

		if (num_of_extbc > 0 && num_extbc_cells > 0)
		{
#ifdef ACTIVE_GPU
			Kernels::compute_extbc_values<<< (num_extbc_cells + THREAD_BLOCK - 1) / THREAD_BLOCK, THREAD_BLOCK, 0, streams >>> (num_extbc_cells, rows, cols, global_dt, device_vec[H], device_vec[QX], device_vec[QY], device_vec[DEM], device_vec[NMAN], device_vec_int[BCRELATIVEINDEX], device_vec_int[BCTYPE], device_vec_int[BCINDEXSTART], device_vec_int[BCNROWSVARS], device_vec[EXTBCV1], device_vec[EXTBCV2], simtime, rank, size);
#else
			Kernels::compute_extbc_values(num_extbc_cells, rows, cols, global_dt, host_vec[H], host_vec[QX], host_vec[QY], host_vec[DEM], host_vec[NMAN], host_vec_int[BCRELATIVEINDEX], host_vec_int[BCTYPE], host_vec_int[BCINDEXSTART], host_vec_int[BCNROWSVARS], host_vec[EXTBCV1], host_vec[EXTBCV2], simtime, rank, size);

#endif

		}

#ifdef ACTIVE_GPU
		Kernels::wet_dry << <(rows*cols + THREAD_BLOCK - 1) / THREAD_BLOCK, THREAD_BLOCK, 0, streams >> > (rows*cols, rows, cols, global_dt, device_vec[H], device_vec[QX], device_vec[QY], device_vec[DEM], device_vec[MAXH], arglist.hextra,size);
#else
		Kernels::wet_dry(rows*cols, rows, cols, global_dt, host_vec[H], host_vec[QX], host_vec[QY], host_vec[DEM], host_vec[MAXH], arglist.hextra,size);
#endif



		if (size > 1)
		{
#ifdef ACTIVE_GPU
			Kernels::halo_copy_from_gpu << <(2 * cols*GHOST_CELL_PADDING + THREAD_BLOCK - 1) / THREAD_BLOCK, THREAD_BLOCK, 0, streams >> > (2 * cols*GHOST_CELL_PADDING, rows, cols, device_vec[H], device_vec[QX], device_vec[QY], device_vec[HALO]);

			if (arglist.gpu_direct_flag)
			{
				cudaStreamSynchronize(streams);
				st.stop(COMPUTE_TIME);

				st.start(BALANCING_MPI_TIME);
				st.start(MPI_TIME);
				MpiUtils::exchange(device_vec[HALO], 12*GHOST_CELL_PADDING, cols, rank, size, USE_HALO);
				st.stop(MPI_TIME);
				st.stop(BALANCING_MPI_TIME);

				st.start(COMPUTE_TIME);
			}
			else
			{
				cudaMemcpyAsync(host_vec[HALO], device_vec[HALO], nbytes_halo, cudaMemcpyDeviceToHost, streams);
				cudaStreamSynchronize(streams);
				st.stop(COMPUTE_TIME);

				st.start(BALANCING_MPI_TIME);
				st.start(MPI_TIME);
				MpiUtils::exchange(host_vec[HALO], 12*GHOST_CELL_PADDING, cols, rank, size, USE_HALO);
				st.stop(MPI_TIME);
				st.stop(BALANCING_MPI_TIME);

				st.start(COMPUTE_TIME);
				cudaMemcpyAsync(device_vec[HALO], host_vec[HALO], nbytes_halo, cudaMemcpyHostToDevice, streams);
			}

			Kernels::halo_copy_to_gpu << <(2 * cols*GHOST_CELL_PADDING + THREAD_BLOCK - 1) / THREAD_BLOCK, THREAD_BLOCK, 0, streams >> > (2 * cols*GHOST_CELL_PADDING, rows, cols, device_vec[H], device_vec[QX], device_vec[QY], device_vec[HALO]);
#else
			Kernels::halo_copy_from_gpu(2 * cols*GHOST_CELL_PADDING, rows, cols, host_vec[H], host_vec[QX], host_vec[QY], host_vec[HALO]);
			st.stop(COMPUTE_TIME);

			st.start(BALANCING_MPI_TIME);
			st.start(MPI_TIME);
			MpiUtils::exchange(host_vec[HALO], 12*GHOST_CELL_PADDING, cols, rank, size, USE_HALO);
			st.stop(MPI_TIME);
			st.stop(BALANCING_MPI_TIME);

			st.start(COMPUTE_TIME);
			Kernels::halo_copy_to_gpu(2 * cols*GHOST_CELL_PADDING, rows, cols, host_vec[H], host_vec[QX], host_vec[QY], host_vec[HALO]);
#endif

#ifdef ACTIVE_GPU
			Kernels::wet_dry_qy_halo << <(2 * cols*GHOST_CELL_PADDING + THREAD_BLOCK - 1) / THREAD_BLOCK, THREAD_BLOCK, 0, streams >> > (2 * cols*GHOST_CELL_PADDING, rows, cols, device_vec[H], device_vec[QY], device_vec[DEM], arglist.hextra);
#else
			Kernels::wet_dry_qy_halo(2 * cols*GHOST_CELL_PADDING, rows, cols, host_vec[H], host_vec[QY], host_vec[DEM], arglist.hextra);
#endif

		}


		st.stop(COMPUTE_TIME);
	}


	template<typename T>
	void triton<T>::new_domain_decomposition()
	{
			
		if(MPI_time_based_domain_decomposition()){;

			reset_arrays();

			process_source_locations();

			process_observation_cells();
			
			process_boundary_condition();

			partition_matrix_files_dynamic();
			
			process_runoff();
			create_host_aux_vectors();
			
			create_host_vectors();
			create_device_vectors();

			//a call to out.init is again neccessary to set the output configuration
			out.init(rows, cols, rank, size, project_dir, arglist.outfile_pattern, arglist.time_series_flag, cfg_content, arglist.output_option);

			if (arglist.time_series_flag)
			{
				out.init_time_series(num_of_obs_points, arglist.observation_x_loc.size(), relative_obs_index, observation_cells, observation_cells_global);

			}
			

		}


	}
	
	template<typename T>
	int triton<T>::MPI_time_based_domain_decomposition()
	{
		int *dyn_rows = new int[size];
		double *mpi_time_all = new double[size];
		double sumMPI;
		double mpi_time = st.get_custom_time(BALANCING_MPI_TIME);
		int sum_rows;
		int flag=0;


		MPI_Gather(&mpi_time, 1, MPI_DATA_TYPE, &mpi_time_all[rank], 1, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);

		if(rank==0){
			sumMPI=0.0;
			for(int j=0;j<size;j++){
				sumMPI+=mpi_time_all[j];
			}
			
			sum_rows=0;
			for(int j=0;j<size;j++){
				double factor = mpi_time_all[j]*size/sumMPI - 1.0;
				if(fabs(factor)>0.05){ //greater than 5%
					flag=1;
				}
				//0.1 is a weight factor to enforce gradual changes in domain size from one step to another
				dyn_rows[j]=pd.part_dims[j].first + (int) floor((factor)*pd.rows/size*0.1);
				//at least 2 real rows per subdomain 
				dyn_rows[j]=max(dyn_rows[j],2+2*GHOST_CELL_PADDING);
				sum_rows+=dyn_rows[j]-2*GHOST_CELL_PADDING;
			}
			
			if(sum_rows<=pd.rows){
				int rem = pd.rows - sum_rows;
				if(rem>0){
					for(int j=0;j<rem;j++){
						dyn_rows[j]++;
					}
				}
			}else{
				std::cerr << ERROR "Re-partitioning algorithm is wrong" << std::endl;
				exit(EXIT_FAILURE);	
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD); 
		
		if(flag==0){
			st.restart(BALANCING_MPI_TIME);
			return 0; 
		}

		MPI_Bcast(dyn_rows, size, MPI_INT, 0, MPI_COMM_WORLD); 


		for(int i=0;i<pd.size;i++){
			pd.part_dims[i].first=dyn_rows[i];
		}

		st.restart(BALANCING_MPI_TIME);
		
		delete[] dyn_rows;
		delete[] mpi_time_all;
		

		return 1;

	}


	template<typename T>
	void triton<T>::reset_arrays()
	{

		delete[] host_vec[EXTBCV2];
		delete[] host_vec[EXTBCV1];
		delete[] host_vec[RUNIN];
		delete[] host_vec[HYGV];
		delete[] host_vec[HYGT];
		delete[] host_vec[DT];
		delete[] host_vec[HALO];
		delete[] host_vec[SQRTH];
		delete[] host_vec[RHSQY1];
		delete[] host_vec[RHSQY0];
		delete[] host_vec[RHSQX1];
		delete[] host_vec[RHSQX0];
		delete[] host_vec[RHSH1];
		delete[] host_vec[RHSH0];

		delete[] host_vec_int[BCNROWSVARS];
		delete[] host_vec_int[BCINDEXSTART];
		delete[] host_vec_int[BCTYPE];
		delete[] host_vec_int[BCRELATIVEINDEX];
		delete[] host_vec_int[RUNID];
		delete[] host_vec_int[SRCP];
		
#ifdef ACTIVE_GPU

		//not neccessary since we are inside output so we already copied this data to the CPU
		/*cudaMemcpyAsync(host_vec[H], device_vec[H], nbytes, cudaMemcpyDeviceToHost, streams);
		cudaMemcpyAsync(host_vec[QX], device_vec[QX], nbytes, cudaMemcpyDeviceToHost, streams);
		cudaMemcpyAsync(host_vec[QY], device_vec[QY], nbytes, cudaMemcpyDeviceToHost, streams);
		if (arglist.max_value_print_option.size() > 0)
		{
			cudaMemcpyAsync(host_vec[MAXH], device_vec[MAXH], nbytes, cudaMemcpyDeviceToHost, streams);
		}
		cudaStreamSynchronize(streams);*/

		cudaStreamDestroy(streams);
		while (!device_vec.empty())
		{
			cudaFree(device_vec.back());
			device_vec.pop_back();
		}
		while (!device_vec_int.empty())
		{
			cudaFree(device_vec_int.back());
			device_vec_int.pop_back();
		}

#endif

	}


	template<typename T>
	void triton<T>::partition_matrix_files_dynamic()
	{		
		if(strcmp(arglist.input_option.c_str(), "SEQ")==0){ //sequential
			sub_dem.resize(1,1);
			sub_nin.resize(1,1);

			sub_dem = MpiUtils::scatter_exchange(dem.get_data(), pd, rank);
			sub_nin = MpiUtils::scatter_exchange(nin.get_data(), pd, rank);
			
			if(arglist.runoff_map.size() > 0)
			{
				sub_rin = MpiUtils::scatter_exchange_int(rin.get_data(), pd, rank);
			}

		}else{
			//gather dem
			if (rank == 0)
			{
				MPI_Gatherv(sub_dem.get_address_at(0, 0), out.cur_proc_data_size, MPI_DATA_TYPE, out.total_data_arr, out.recvcounts, out.displs, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);
			}
			else
			{
				MPI_Gatherv(sub_dem.get_address_at(GHOST_CELL_PADDING, 0), out.cur_proc_data_size, MPI_DATA_TYPE, out.total_data_arr, out.recvcounts, out.displs, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);
			}
			sub_dem.resize(1,1);
			sub_dem = MpiUtils::scatter_exchange(out.total_data_arr, pd, rank);
			
			//gather nin
			if (rank == 0)
			{
				MPI_Gatherv(sub_nin.get_address_at(0, 0), out.cur_proc_data_size, MPI_DATA_TYPE, out.total_data_arr, out.recvcounts, out.displs, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);
			}
			else
			{
				MPI_Gatherv(sub_nin.get_address_at(GHOST_CELL_PADDING, 0), out.cur_proc_data_size, MPI_DATA_TYPE, out.total_data_arr, out.recvcounts, out.displs, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);
			}
			sub_nin.resize(1,1);
			sub_nin = MpiUtils::scatter_exchange(out.total_data_arr, pd, rank);

			if(arglist.runoff_map.size() > 0)
			{	
				//gather rmap
				if (rank == 0)
				{
					MPI_Gatherv(sub_rin.get_address_at(0, 0), out.cur_proc_data_size, MPI_INTEGER, out.total_data_arr_int, out.recvcounts, out.displs, MPI_INTEGER, 0, MPI_COMM_WORLD);
				}
				else
				{
					MPI_Gatherv(sub_rin.get_address_at(GHOST_CELL_PADDING, 0), out.cur_proc_data_size, MPI_INTEGER, out.total_data_arr_int, out.recvcounts, out.displs, MPI_INTEGER, 0, MPI_COMM_WORLD);
				}
				sub_rin.resize(1,1);
				sub_rin = MpiUtils::scatter_exchange_int(out.total_data_arr_int, pd, rank);
			}


		}
		
		rows = sub_dem.get_num_rows();
		cols = sub_dem.get_num_cols();
			
		sub_dem.set_nrows(sub_dem.get_num_rows());
		sub_dem.set_ncols(sub_dem.get_num_cols());
		sub_dem.set_cell_size(dem.get_cell_size());
		sub_dem.set_xll_corner(dem.get_xll_corner());
		sub_dem.set_yll_corner(dem.get_yll_corner());
		sub_dem.set_no_data_value(dem.get_no_data_value());


		//gather H
		if (rank == 0)
		{
			MPI_Gatherv(sub_hin.get_address_at(0, 0), out.cur_proc_data_size, MPI_DATA_TYPE, out.total_data_arr, out.recvcounts, out.displs, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);
		}
		else
		{
			MPI_Gatherv(sub_hin.get_address_at(GHOST_CELL_PADDING, 0), out.cur_proc_data_size, MPI_DATA_TYPE, out.total_data_arr, out.recvcounts, out.displs, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);
		}

		sub_hin.resize(1,1);
		sub_hin = MpiUtils::scatter_exchange(out.total_data_arr, pd, rank);


		//gather QX
		if (rank == 0)
		{
			MPI_Gatherv(sub_qxin.get_address_at(0, 0), out.cur_proc_data_size, MPI_DATA_TYPE, out.total_data_arr, out.recvcounts, out.displs, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);
		}
		else
		{
			MPI_Gatherv(sub_qxin.get_address_at(GHOST_CELL_PADDING, 0), out.cur_proc_data_size, MPI_DATA_TYPE, out.total_data_arr, out.recvcounts, out.displs, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);
		}
		sub_qxin.resize(1,1);
		sub_qxin = MpiUtils::scatter_exchange(out.total_data_arr, pd, rank);
		
		//gather QY
		if (rank == 0)
		{
			MPI_Gatherv(sub_qyin.get_address_at(0, 0), out.cur_proc_data_size, MPI_DATA_TYPE, out.total_data_arr, out.recvcounts, out.displs, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);
		}
		else
		{
			MPI_Gatherv(sub_qyin.get_address_at(GHOST_CELL_PADDING, 0), out.cur_proc_data_size, MPI_DATA_TYPE, out.total_data_arr, out.recvcounts, out.displs, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);
		}
		sub_qyin.resize(1,1);
		sub_qyin = MpiUtils::scatter_exchange(out.total_data_arr, pd, rank);


		//gather MAXH
		if (rank == 0)
		{
			MPI_Gatherv(sub_max_value_h.get_address_at(0, 0), out.cur_proc_data_size, MPI_DATA_TYPE, out.total_data_arr, out.recvcounts, out.displs, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);
		}
		else
		{
			MPI_Gatherv(sub_max_value_h.get_address_at(GHOST_CELL_PADDING, 0), out.cur_proc_data_size, MPI_DATA_TYPE, out.total_data_arr, out.recvcounts, out.displs, MPI_DATA_TYPE, 0, MPI_COMM_WORLD);
		}
		sub_max_value_h.resize(1,1);
		sub_max_value_h = MpiUtils::scatter_exchange(out.total_data_arr, pd, rank);

		if(rank == 0){
			std::cerr << OK "Data has been re-partitioned" << std::endl;
		}

	}


}

#endif
