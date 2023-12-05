/** @file config_utils.h
 *  @brief Header containing the ConfigUtils class
 *
 *  This contains the subroutines and eventually any 
 *  macros, constants, etc. needed for ConfigUtils class
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

#ifndef CONFIG_UTILS_H
#define CONFIG_UTILS_H

#include "string_utils.h"

namespace ConfigUtils
{
	template<typename T>
	struct arguments	/**< Structure to contain all arguments extracted from configuration (cfg) file. */
	{
		bool
		time_increment_fixed,	/**< Flag to indicate time step size characteristics. True = Constant time step size, False = Variable time step size. */
		time_series_flag,	/**< Flag to allow time series output. True = Output time series, False = Don't output time series. */
		gpu_direct_flag,	/**< Flag to allow GPU-Direct use. True = Use GPU-Direct, False = Don't use GPU-Direct. */
		open_boundaries; /**<Flag to impose open or closed boundaries. By default, open_boundaries=0*/ 

		int
		checkpoint_id,	/**< Use for hot start. If 0 then that means a clean start. Greater than 0 means start from that specific checkpoint. */
		num_sources,	/**< The total number of flow locations in Hygrograph. If there are no flow locations then 0 is allowed. */
		num_runoffs,	/**< The total number of Runoffs. */
		num_extbc,	/**< The total number of External boundary cells group. Each group can contain one or multiple cells. */
		it_count,	/**< The total number of iterations up to a specific point. 0 in case of a clean start, greater than 0 otherwise. */
		factor_interval_domain_decomposition;	/**< Factor applied to the print interval time to check for domain decomposition. */

		T
		time_step,	/**< Indicates the time step size. Time step size determines the time for the next computation. */
		sim_start_time,	/**< Starting time point of a simulation. Usually 0 for a new simulation. */
		sim_duration,	/**< Finishing time point of a simulation. Regardless of the starting point, simulation always ends at this point. */
		print_interval,	/**< Use for outputting files. After every defined print interval time, the program will save outputs in an external file. */
		courant,	/**< Represents Courant number. */
		const_mann,	/**< Constant manning value to use in every cell in case of no external manning file is provided. */
		hextra;	/**< Represents a the minimum water depth tolerance */

		std::string
		outfile_pattern,	/**< Output file directory and name pattern. */
		hydrograph_filename,	/**< Directory of the Hygrograph file to use. */
		runoff_filename,	/**< Directory of the Runoff file to use. */
		print_option,	/**< Use to determine output types. h to output just the h (depth), huv to output all h (depth),u and v (velocities). */
		max_value_print_option,	/**< Use to determine maximum value of each cells output types. h to output just the h (depth). */
		input_format,	/**< Expected input file format. BIN for binary file or ASC for ascii file. */
		input_option,	/**< Strategy to use for input files. PAR for parallel input or SEQ for sequential inp. PAR reads each MPI partitions subdomain in separate files and SEQ reads the whole domain from one file. Applied to all raster formats. By default, SEQ is considered.*/
		output_format,	/**< Expected output file format. BIN for binary file or ASC for ascii file. */
		output_option,	/**< Strategy to use for outputting into files. PAR for parallel outputs or SEQ for sequential outputs. PAR saves each MPI partitions subdomain in separate files and SEQ saves the whole domain into one file. */
		dem_filename,	/**< Directory of the DEM file to use. */
		header_filename,	/**< Directory of the header file to use (in case of parallel reading). */
		src_loc_file,	/**< Directory of the file that contains the information of all flow locations. */
		runoff_map,	/**< Directory of the Runoff map to use. */
		observation_loc_file,	/**< Directory of the file that contains the information of all cells to observe and generate time series output. */
		extbc_file,	/**< Directory of the External boundary condition file to use. */
		extbc_dir,	/**< Parent directory of the External boundary condition files. */
		h_infile,	/**< Initial water depth file directory. */
		qx_infile,	/**< Initial flux in x direction file directory. */
		qy_infile,	/**< Initial flux in y direction file directory. */
		n_infile,	/**< Directory of the manning file to use. */
		domain_decomposition;	/**< Domain decomposition. Options are static or dynamic. Static by default*/


		std::vector<T>
		src_x_loc,	/**< Vector to hold all the Longitude value of all the flow locations serially. */
		src_y_loc,	/**< Vector to hold all the Latitude value of all the flow locations serially. */
		observation_x_loc,	/**< Vector to hold all the Longitude value of all the observation cells. */
		observation_y_loc,	/**< Vector to hold all the Latitude value of all the observation cells. */
		extbc_x1_loc,	/**< Vector to hold all the Longitude value of the starting cell of an external boundary condition. */
		extbc_y1_loc,	/**< Vector to hold all the Latitude value of the starting cell of an external boundary condition. */
		extbc_x2_loc,	/**< Vector to hold all the Longitude value of the ending cell of an external boundary condition. */
		extbc_y2_loc;	/**< Vector to hold all the Latitude value of the ending cell of an external boundary condition. */

		std::vector<int> extbc_bctype;	/**< Contains all external boundary condition type serially. */
		std::vector<std::string> extbc_fname;	/**< Contains all external boundary condition file name serially. */

	};


/** @brief It calculates the corresponding value of each attribute name from the contents of the configuration (cfg) file.
*
*  @param x attibute name
*  @param y contents of cfg file
*  @param d default value
*  @return The corresponding value
*/
	std::string argsd(std::string x, std::map<std::string, std::string> y, std::string d);


/** @brief It calculates the corresponding value of each attribute name from the contents of the configuration (cfg) file without any default value.
*
*  @param x attibute name
*  @param y contents of cfg file
*  @return The corresponding value
*/
	std::string args(std::string x, std::map<std::string, std::string> y);


/** @brief It extracts the whole configuration string and constructs an attribute key-value mapping.
*
*  @param cfg_content cfg file content
*  @return Attribute key value mapping
*/
	std::map<std::string, std::string> parse_cfg(std::string cfg_content);


/** @brief It extracts each flow location and observation cells Longitude and Latitude value and constructs a (x,y) location mapping.
*
*  @param filename file to parse
*  @param type determine flow location of observation
*  @return (x,y) location mapping
*/
	std::map<std::string, std::string> parse_src_location(std::string filename, int type);


/** @brief It extracts each external boundary condition file and constructs an attribute key-value mapping.
*
*  @param filename file to parse
*  @param dir parent directory of filename
*  @return Attribute key value mapping
*/
	std::map<std::string, std::string> parse_extbc_file(std::string filename, std::string dir);


/** @brief It calculates all argument values and constructs struct arguments object.
*
*  @param cfg file to parse
*  @return The arguments object contating all argument
*/
	template<typename T>
	arguments<T> get_args(std::string cfg);


/** @brief It reads a configuration (cfg) file and constructs a string of the whole file.
*
*  @param filepath file to read
*  @return Contents as a string
*/
	std::string file_content_to_string(std::string filepath);


/** @brief It computes the root directory from the full path, shortening it out when a backslash is found.
*
*  @param path The full path
*  @return A string with the project directory
*/
	std::string get_root_dir(const char* path);

/** @brief It reads the number of rows from the output file when checkpoint is enabled.
*
*  @param project_dir String containing the project directory
*  @param dyn_rows Array of size "number of ranks" that will contain the number of rows
*  @param checkpoint_id Checkpoint id
*/
	void read_and_parse_checkpoint_partition(std::string project_dir, int *dyn_rows, int checkpoint_id);


	std::string argsd(std::string x, std::map<std::string, std::string> y, std::string d)
	{
		if(y.find(x) != y.end())
		{
			return (y.find(x))->second;
		}
		else
		{
			return d;
		}
	}


	std::string args(std::string x, std::map<std::string, std::string> y)
	{
		return argsd(x, y, "");
	}


	std::map<std::string, std::string> parse_cfg(std::string cfg_content)
	{
		std::map<std::string, std::string> arglist;
		std::istringstream ifs(cfg_content.c_str());
		std::string line;

		if (!ifs.good())
		{
			std::cerr << "Error parsing configuration content." << std::endl;
			exit(EXIT_FAILURE);
		}

		while (std::getline(ifs, line))
		{
			line = StringUtils::trim(line);

			if (line.size() > 0 && line[0] != '#')
			{
				size_t incomm = line.find('#');
				if ((incomm != std::string::npos) && (line[incomm - 1] != '\'') && (line[incomm - 1] != '\"'))
				{
					line.erase(line.begin() + incomm);
				}
				std::vector<std::string> kv = StringUtils::split(line, '=');
				std::string key = kv[0];
				std::string value = kv[1];

				if (value[0] == '"' && value[value.size() - 1] == '"')
				{
					value.erase(value.begin());
					value.erase(value.end() - 1);
				}

				if (kv.size() > 2)
				{
					std::vector<std::string>::iterator it = kv.begin() + 2;
					for (; it != kv.end(); it++)
					{
						value += (*it);
					}
				}
				arglist.insert(std::pair<std::string, std::string>(key, value));
			}
		}

		return arglist;
	}


	std::map<std::string, std::string> parse_src_location(std::string filename, int type)
	{
		std::map<std::string, std::string> arglist;
		std::ifstream ifs(filename.c_str());
		std::string line;

		if (!ifs.good())
		{
			std::cerr << ERROR "Error reading file: " << filename << std::endl;
			exit(EXIT_FAILURE);
		}
		std::string src_x_loc = "", src_y_loc = "";

		while (std::getline(ifs, line))
		{
			line = StringUtils::trim(line);

			if (line.size() > 0 && line[0] != '%')
			{
				size_t incomm = line.find('#');
				if ((incomm != std::string::npos) && (line[incomm - 1] != '\'') && (line[incomm - 1] != '\"'))
				{
					line.erase(line.begin() + incomm);
				}
				std::vector<std::string> kv = StringUtils::split(line, ',');
				std::string x_value = kv[0];
				std::string y_value = kv[1];

				if (x_value[0] == '"' && x_value[x_value.size() - 1] == '"')
				{
					x_value.erase(x_value.begin());
					x_value.erase(x_value.end() - 1);
				}
				if (y_value[0] == '"' && y_value[y_value.size() - 1] == '"')
				{
					y_value.erase(y_value.begin());
					y_value.erase(y_value.end() - 1);
				}

				src_x_loc = src_x_loc + "," + x_value;
				src_y_loc = src_y_loc + "," + y_value;
			}
		}

		src_x_loc = src_x_loc.substr(1);
		src_y_loc = src_y_loc.substr(1);

		if (type == SRC_LOCATION)
		{
			arglist.insert(std::pair<std::string, std::string>("src_x_loc", src_x_loc));
			arglist.insert(std::pair<std::string, std::string>("src_y_loc", src_y_loc));
		}
		else if (type == OBSERVATION_LOCATION)
		{
			arglist.insert(std::pair<std::string, std::string>("observation_x_loc", src_x_loc));
			arglist.insert(std::pair<std::string, std::string>("observation_y_loc", src_y_loc));
		}

		ifs.close();

		return arglist;
	}


	std::map<std::string, std::string> parse_extbc_file(std::string filename, std::string dir)
	{
		std::map<std::string, std::string> arglist;
		std::ifstream ifs(filename.c_str());
		std::string line;

		if (!ifs.good())
		{
			std::cerr << ERROR "Error reading file: " << filename << std::endl;
			exit(EXIT_FAILURE);
		}

		std::string bctype = "";
		std::string src_x1_loc = "", src_y1_loc = "";
		std::string src_x2_loc = "", src_y2_loc = "";
		std::string bcfname = "";

		while (std::getline(ifs, line))
		{
			line = StringUtils::trim(line);

			if (line.size() > 0 && line[0] != '%')
			{
				size_t incomm = line.find('#');
				if ((incomm != std::string::npos) && (line[incomm - 1] != '\'') && (line[incomm - 1] != '\"'))
				{
					line.erase(line.begin() + incomm);
				}
				std::vector<std::string> kv = StringUtils::split(line, ',');
				std::string bctype_value = kv[0];
				std::string x1_value = kv[1];
				std::string y1_value = kv[2];
				std::string x2_value = kv[3];
				std::string y2_value = kv[4];

				int auxtype=std::stoi(kv[0]);
				std::string bcfname_value ="";
				if(auxtype!=0)
				{
					bcfname_value = kv[5];
				}
				else
				{
					bcfname_value = "0.0";
				}

				if (bctype_value[0] == '"' && bctype_value[bctype_value.size() - 1] == '"')
				{
					bctype_value.erase(bctype_value.begin());
					bctype_value.erase(bctype_value.end() - 1);
				}

				if (x1_value[0] == '"' && x1_value[x1_value.size() - 1] == '"')
				{
					x1_value.erase(x1_value.begin());
					x1_value.erase(x1_value.end() - 1);
				}
				if (y1_value[0] == '"' && y1_value[y1_value.size() - 1] == '"')
				{
					y1_value.erase(y1_value.begin());
					y1_value.erase(y1_value.end() - 1);
				}

				if (x2_value[0] == '"' && x2_value[x2_value.size() - 1] == '"')
				{
					x2_value.erase(x2_value.begin());
					x2_value.erase(x2_value.end() - 1);
				}
				if (y2_value[0] == '"' && y2_value[y2_value.size() - 1] == '"')
				{
					y2_value.erase(y2_value.begin());
					y2_value.erase(y2_value.end() - 1);
				}

				if (bcfname_value[0] == '"' && bcfname_value[bcfname_value.size() - 1] == '"')
				{
					bcfname_value.erase(bcfname_value.begin());
					bcfname_value.erase(bcfname_value.end() - 1);
				}

				bctype = bctype + "," + bctype_value;
				src_x1_loc = src_x1_loc + "," + x1_value;
				src_y1_loc = src_y1_loc + "," + y1_value;
				src_x2_loc = src_x2_loc + "," + x2_value;
				src_y2_loc = src_y2_loc + "," + y2_value;

				if(auxtype==1)
				{
					bcfname = bcfname + "," + dir + "/" + bcfname_value;
				}
				else   //case 0, 2, 3
				{
					bcfname = bcfname + "," + bcfname_value;
				}
			}
		}

		bctype = bctype.substr(1);
		src_x1_loc = src_x1_loc.substr(1);
		src_y1_loc = src_y1_loc.substr(1);
		src_x2_loc = src_x2_loc.substr(1);
		src_y2_loc = src_y2_loc.substr(1);
		bcfname = bcfname.substr(1);

		arglist.insert(std::pair<std::string, std::string>("extbc_bctype", bctype));
		arglist.insert(std::pair<std::string, std::string>("extbc_x1_loc", src_x1_loc));
		arglist.insert(std::pair<std::string, std::string>("extbc_y1_loc", src_y1_loc));
		arglist.insert(std::pair<std::string, std::string>("extbc_x2_loc", src_x2_loc));
		arglist.insert(std::pair<std::string, std::string>("extbc_y2_loc", src_y2_loc));
		arglist.insert(std::pair<std::string, std::string>("extbc_fname", bcfname));

		ifs.close();

		return arglist;
	}


	template<typename T>
	arguments<T> get_args(std::string cfg)
	{
		std::map<std::string, std::string> argmap = parse_cfg(cfg);

		arguments<T> arglist;

		arglist.open_boundaries=0; //by default
		arglist.input_option = "SEQ"; //by default

		arglist.outfile_pattern = argsd("outfile_pattern", argmap, "");
		arglist.hydrograph_filename = args("hydrograph_filename", argmap);
		arglist.runoff_filename = args("runoff_filename", argmap);
		arglist.dem_filename = args("dem_filename", argmap); //in the case of PAR input, the dem file should include only the name (without extension and _)
		arglist.src_loc_file = argsd("src_loc_file", argmap, "");
		arglist.runoff_map = argsd("runoff_map", argmap, ""); //in the case of PAR input, the rmap file should include only the name (without extension and _)
		arglist.observation_loc_file = argsd("observation_loc_file", argmap, "");
		arglist.extbc_file = argsd("extbc_file", argmap, "");
		arglist.extbc_dir = argsd("extbc_dir", argmap, "");

		arglist.const_mann = atof((args("const_mann", argmap)).c_str());
		arglist.time_step = atof((args("time_step", argmap)).c_str());
		arglist.time_increment_fixed = atoi((args("time_increment_fixed", argmap)).c_str());
		arglist.time_series_flag = atoi((args("time_series_flag", argmap)).c_str());
		arglist.courant = atof((args("courant", argmap)).c_str());
		arglist.hextra = atof((args("hextra", argmap)).c_str());
		arglist.gpu_direct_flag = atoi((args("gpu_direct_flag", argmap)).c_str());
		arglist.open_boundaries = atoi((args("open_boundaries", argmap)).c_str());

		arglist.num_sources = atoi((args("num_sources", argmap)).c_str());
		arglist.num_runoffs = atoi((args("num_runoffs", argmap)).c_str());
		arglist.num_extbc = atoi((args("num_extbc", argmap)).c_str());
		arglist.checkpoint_id = atoi((argsd("checkpoint_id", argmap, "0")).c_str());
		arglist.it_count = atoi((argsd("it_count", argmap, "0")).c_str());
		arglist.print_option = StringUtils::tolower(args("print_option", argmap));
		arglist.max_value_print_option = argsd("max_value_print_option", argmap, "");
		arglist.input_format = args("input_format", argmap);
		arglist.input_option = args("input_option", argmap);
		arglist.output_format = args("output_format", argmap);
		arglist.output_option = args("output_option", argmap);

		arglist.h_infile = argsd("h_infile", argmap, "");
		arglist.qx_infile = argsd("qx_infile", argmap, "");
		arglist.qy_infile = argsd("qy_infile", argmap, "");
		arglist.n_infile = argsd("n_infile", argmap, "");  //in the case of PAR input, the n_infile file should include only the name (without extension and _)

		arglist.domain_decomposition = args("domain_decomposition", argmap);
		arglist.factor_interval_domain_decomposition = atoi((args("factor_interval_domain_decomposition", argmap)).c_str());

		arglist.sim_start_time = atof((args("sim_start_time", argmap)).c_str());
		arglist.sim_duration = atof((args("sim_duration", argmap)).c_str());
		arglist.print_interval = atof((args("print_interval", argmap)).c_str());

		if(arglist.num_sources > 0)
		{
			std::map<std::string, std::string> src_map = parse_src_location(arglist.src_loc_file, SRC_LOCATION);
			arglist.src_x_loc = StringUtils::vecstr_to_vecflt<T>(StringUtils::split((args("src_x_loc", src_map)), ','));
			arglist.src_y_loc = StringUtils::vecstr_to_vecflt<T>(StringUtils::split((args("src_y_loc", src_map)), ','));
		}

		if (arglist.time_series_flag)
		{
			std::map<std::string, std::string> observation_loc_map = parse_src_location(arglist.observation_loc_file, OBSERVATION_LOCATION);
			arglist.observation_x_loc = StringUtils::vecstr_to_vecflt<T>(StringUtils::split((args("observation_x_loc", observation_loc_map)), ','));
			arglist.observation_y_loc = StringUtils::vecstr_to_vecflt<T>(StringUtils::split((args("observation_y_loc", observation_loc_map)), ','));
		}

		//external BC
		if (arglist.num_extbc > 0)
		{
			std::map<std::string, std::string> extbc_map = parse_extbc_file(arglist.extbc_file,arglist.extbc_dir);
			arglist.extbc_bctype = StringUtils::vecstr_to_vecint(StringUtils::split((args("extbc_bctype", extbc_map)), ','));
			arglist.extbc_x1_loc = StringUtils::vecstr_to_vecflt<T>(StringUtils::split((args("extbc_x1_loc", extbc_map)), ','));
			arglist.extbc_y1_loc = StringUtils::vecstr_to_vecflt<T>(StringUtils::split((args("extbc_y1_loc", extbc_map)), ','));
			arglist.extbc_x2_loc = StringUtils::vecstr_to_vecflt<T>(StringUtils::split((args("extbc_x2_loc", extbc_map)), ','));
			arglist.extbc_y2_loc = StringUtils::vecstr_to_vecflt<T>(StringUtils::split((args("extbc_y2_loc", extbc_map)), ','));
			arglist.extbc_fname = StringUtils::split((args("extbc_fname", extbc_map)), ',');
		}
		
		if(strcmp(arglist.input_option.c_str(), "")==0){ //no input_option provided (for legacy version)
			arglist.input_option = "SEQ"; //by default
		}

		if(strcmp(arglist.input_option.c_str(), "PAR")==0){
			arglist.header_filename = args("header_filename", argmap);
		}

		return arglist;
	}


	std::string file_content_to_string(std::string filepath)
	{
		std::string file_content;
		std::ifstream input(filepath);
		if (input.is_open())
		{
			input.seekg(0, std::ios::end);
			file_content.reserve(input.tellg());
			input.seekg(0, std::ios::beg);
			file_content.assign((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
			input.close();
		}
		else
		{
			std::cerr << ERROR "Error reading file: " << filepath << std::endl;
			exit(EXIT_FAILURE);
		}

		return file_content;
	}


	std::string get_root_dir(const char* path)
	{
		std::string spath(path);
		int found = spath.find_last_of("/\\");
		spath = spath.substr(0, found);

		return spath.substr(0, (spath.find_last_of("/\\")));
	}

	void read_and_parse_checkpoint_partition(std::string project_dir, int *dyn_rows, int checkpoint_id)
	{
		
		std::string outdir = project_dir + "/" + OUTPUT_DIR + "/domain_decomposition/";
		std::string filename = outdir + "domain_decomposition" + std::to_string(checkpoint_id) + ".txt";
		
		std::ifstream ifs(filename.c_str());
		std::string line;

		if (!ifs.good())
		{
			std::cerr << ERROR "Error reading file: " << filename << std::endl;
			exit(EXIT_FAILURE);
		}

		while (std::getline(ifs, line))
		{
			line = StringUtils::trim(line);

			if (line.size() > 0 && line[0] != '%')
			{
				std::vector<std::string> kv = StringUtils::split(line, ',');
				std::string rank_id = kv[0];
				std::string nrows = kv[1];
				std::string ncols = kv[2];
				int rank_int=std::stoi(kv[0]);
				dyn_rows[rank_int]=std::stoi(kv[1]);
			}
		}

		ifs.close();

	}

}

#endif
