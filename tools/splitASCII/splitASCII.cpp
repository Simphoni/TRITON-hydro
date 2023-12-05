#include <iostream>
#include <fstream>
#include <limits> 
#include <sstream>
#include <vector>
#include <dirent.h>
#include <sys/stat.h>
#include <mpi.h>
#include <math.h>
#include <algorithm> // Add this line to include the <algorithm> header


const int NFILES = 7;
const std::string IS_MANN = "YES";
const std::string IS_RMAP = "YES";

const std::string TRITON_DIR = "/home/mario/fork_tritonmpi/tritonmpi";
const std::string INPUT_DEM = TRITON_DIR + "/input/dem/asc/case03.dem";
const std::string INPUT_MANN = TRITON_DIR + "/input/mann/asc/case03.mann";
const std::string INPUT_RMAP = TRITON_DIR + "/input/runoff/case03_runoff.rmap";


std::vector<std::string>& split(const std::string& s, char delim, std::vector<std::string>& elems)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim))
    {
        elems.emplace_back(std::move(item));
    }
    return elems;
}

std::vector<std::string> split(const std::string& s, char delim)
{
    std::vector<std::string> elems;
    elems.reserve(std::count(s.begin(), s.end(), delim) + 1);
    return split(s, delim, elems);
}


void split_dem_to_bin(const std::string& casename_dem, const int start_idx, const int end_idx, const std::vector<int> line_numbers, const long ncols) {
    
	std::ifstream input(INPUT_DEM);
   if (!input.is_open()) {
       std::cerr << "Error opening input DEM file." << std::endl;
       return;
   }

   std::string line;
   int line_count = -6;
	int start_index= ((start_idx == 0) ? 0 : line_numbers[start_idx - 1]);
	int end_index=line_numbers[end_idx - 1];
	int nrows_total=end_index-start_index;

   double *arr = new double [nrows_total*ncols];

   // Skip lines until reaching start_index
   while (line_count < start_index  && std::getline(input, line)) {
   	line_count++;
   }

	long i = 0;
   while (line_count < end_index && std::getline(input, line)) {
   	std::vector<std::string> row = split(line, ' ');
   	for (long j = 0; j < row.size(); j++) {
   		arr[(ncols * i) + j] = (double)atof(row[j].c_str());
    	}
		i++;
		line_count++;
	}

	input.close();

	int sum=0;
	for (int idx = start_idx; idx < end_idx; idx++) {
		start_index = (idx == 0) ? 0 : line_numbers[idx - 1] ;
		end_index = line_numbers[idx];
	 	int nrows_local=end_index-start_index;

		 std::string outfile = casename_dem + "_" + (idx < 10 ? "0" : "") + std::to_string(idx) + ".dem";
		 std::ofstream output(outfile, std::ios::binary);

		 if (!output.is_open()) {
			  std::cerr << "Error opening file: " << outfile << std::endl;
			  return;
		 }
		 double put_rows_value = (double)(nrows_local);
		 double put_cols_value = (double)(ncols);
				
		 output.write((char*) &put_rows_value, sizeof(double));
		 output.write((char*) &put_cols_value, sizeof(double));

		 output.write((char*)&arr[sum], nrows_local*ncols * sizeof(double));
		 output.close();

		 sum+=nrows_local*ncols;
		 std::cout << "Split ASCII DEM file and converted to BIN " << outfile << std::endl;

	 }
	 delete[] arr; 

}


void split_mann_to_bin(const std::string& casename_mann, const int start_idx, const int end_idx, const std::vector<int> line_numbers, const long ncols) {
    
	std::ifstream input(INPUT_MANN);
   if (!input.is_open()) {
       std::cerr << "Error opening input MANN file." << std::endl;
       return;
   }

   std::string line;
   int line_count = 0;
	int start_index= ((start_idx == 0) ? 0 : line_numbers[start_idx - 1]);
	int end_index=line_numbers[end_idx - 1];
	int nrows_total=end_index-start_index;

   double *arr = new double [nrows_total*ncols];

   // Skip lines until reaching start_index
   while (line_count < start_index  && std::getline(input, line)) {
   	line_count++;
   }

	long i = 0;
   while (line_count < end_index && std::getline(input, line)) {
   	std::vector<std::string> row = split(line, ' ');
   	for (long j = 0; j < row.size(); j++) {
   		arr[(ncols * i) + j] = (double)atof(row[j].c_str());
    	}
		i++;
		line_count++;
	}

	input.close();

	int sum=0;
	for (int idx = start_idx; idx < end_idx; idx++) {
		start_index = (idx == 0) ? 0 : line_numbers[idx - 1] ;
		end_index = line_numbers[idx];
	 	int nrows_local=end_index-start_index;

		 std::string outfile = casename_mann + "_" + (idx < 10 ? "0" : "") + std::to_string(idx) + ".mann";
		 std::ofstream output(outfile, std::ios::binary);

		 if (!output.is_open()) {
			  std::cerr << "Error opening file: " << outfile << std::endl;
			  return;
		 }
		 double put_rows_value = (double)(nrows_local);
		 double put_cols_value = (double)(ncols);
				
		 output.write((char*) &put_rows_value, sizeof(double));
		 output.write((char*) &put_cols_value, sizeof(double));

		 output.write((char*)&arr[sum], nrows_local*ncols * sizeof(double));
		 output.close();

		 sum+=nrows_local*ncols;
		 std::cout << "Split ASCII MANN file and converted to BIN " << outfile << std::endl;

	 }
	 delete[] arr; 

}

void split_rmap_to_bin(const std::string& casename_rmap, const int start_idx, const int end_idx, const std::vector<int> line_numbers, const long ncols) {
    
	std::ifstream input(INPUT_RMAP);
   if (!input.is_open()) {
       std::cerr << "Error opening input RMAP file." << std::endl;
       return;
   }

   std::string line;
   int line_count = 0;
	int start_index= ((start_idx == 0) ? 0 : line_numbers[start_idx - 1]);
	int end_index=line_numbers[end_idx - 1];
	int nrows_total=end_index-start_index;

   int *arr = new int [nrows_total*ncols];


   // Skip lines until reaching start_index
   while (line_count < start_index  && std::getline(input, line)) {
   	line_count++;
   }

	long i = 0;
   while (line_count < end_index && std::getline(input, line)) {
   	std::vector<std::string> row = split(line, ' ');
   	for (long j = 0; j < row.size(); j++) {
   		arr[(ncols * i) + j] = atoi(row[j].c_str());
    	}
		i++;
		line_count++;
	}

	input.close();

	int sum=0;
	for (int idx = start_idx; idx < end_idx; idx++) {
		start_index = (idx == 0) ? 0 : line_numbers[idx - 1] ;
		end_index = line_numbers[idx];
	 	int nrows_local=end_index-start_index;

		 std::string outfile = casename_rmap + "_" + (idx < 10 ? "0" : "") + std::to_string(idx) + ".rmap";
		 std::ofstream output(outfile, std::ios::binary);

		 if (!output.is_open()) {
			  std::cerr << "Error opening file: " << outfile << std::endl;
			  return;
		 }
		 int put_rows_value = (int)(nrows_local);
		 int put_cols_value = (int)(ncols);
				
		 output.write((char*) &put_rows_value, sizeof(int));
		 output.write((char*) &put_cols_value, sizeof(int));

		 output.write((char*)&arr[sum], nrows_local*ncols * sizeof(int));
		 output.close();

		 sum+=nrows_local*ncols;
		 std::cout << "Split ASCII RMAP file and converted to BIN " << outfile << std::endl;

	 }
	 delete[] arr; 

}


int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
	 std::string filename;
	 std::string output_dir, parent_dir;
	 std::string filenameWithoutExtension;
	 std::size_t lastSlashPos, lastDotPos;

	 std::ifstream input(INPUT_DEM);
    std::string line;
	 
	 long ncols;
	 long nrows;
	 DIR* dir_;


    if (input.is_open()) {
        std::getline(input, line);
        std::istringstream iss1(line); // Create an input string stream from the line
        iss1.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
        iss1 >> ncols; // Read the value of ncols

        std::getline(input, line);
        std::istringstream iss2(line); // Create an input string stream from the line
        iss2.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
        iss2 >> nrows; // Read the value of nrows

        input.close();
		  if(rank==0){
        	std::cout << "Number of columns: " << ncols << std::endl;
        	std::cout << "Number of rows: " << nrows << std::endl;
    	  }
	 } else {
        std::cerr << "Error opening DEM file." << std::endl;
        MPI_Finalize();
        return 1;
    }


    std::string input_dir_dem = INPUT_DEM.substr(0, INPUT_DEM.find_last_of("/"));
    std::string input_dir_mann = INPUT_MANN.substr(0, INPUT_MANN.find_last_of("/"));
    std::string input_dir_rmap = INPUT_RMAP.substr(0, INPUT_RMAP.find_last_of("/"));


    // DEM: Extract the filename from the path
    lastSlashPos = INPUT_DEM.find_last_of('/');
    filename= INPUT_DEM.substr(lastSlashPos + 1);
    // DEM: Remove the file extension
    lastDotPos = filename.find_last_of('.');
    filenameWithoutExtension = filename.substr(0, lastDotPos);
	 output_dir= input_dir_dem + "/par-bin/";
	 parent_dir = output_dir.substr(0, output_dir.find_last_of("/"));
	 dir_ = opendir(parent_dir.empty() ? "." : parent_dir.c_str());
	 if( !dir_) {
		mkdir(parent_dir.c_str(), S_IRWXU);
	 }
	 else closedir(dir_);
	 std::string casename_dem = output_dir + filenameWithoutExtension;


    // MANN: Extract the filename from the path
    lastSlashPos = INPUT_MANN.find_last_of('/');
    filename= INPUT_MANN.substr(lastSlashPos + 1);
    // MANN: Remove the file extension
    lastDotPos = filename.find_last_of('.');
    filenameWithoutExtension = filename.substr(0, lastDotPos);
	 output_dir= input_dir_mann + "/par-bin/";
	 parent_dir = output_dir.substr(0, output_dir.find_last_of("/"));
	 dir_ = opendir(parent_dir.empty() ? "." : parent_dir.c_str());
	 if( !dir_) {
		mkdir(parent_dir.c_str(), S_IRWXU);
	 }
	 else closedir(dir_);
	 std::string casename_mann = output_dir + filenameWithoutExtension;


    // RMAP: Extract the filename from the path
    lastSlashPos = INPUT_RMAP.find_last_of('/');
    filename= INPUT_RMAP.substr(lastSlashPos + 1);
    // RMAP: Remove the file extension
    lastDotPos = filename.find_last_of('.');
    filenameWithoutExtension = filename.substr(0, lastDotPos);
	 output_dir= input_dir_rmap + "/par-bin/";
	 parent_dir = output_dir.substr(0, output_dir.find_last_of("/"));
	 dir_ = opendir(parent_dir.empty() ? "." : parent_dir.c_str());
	 if( !dir_) {
		mkdir(parent_dir.c_str(), S_IRWXU);
	 }
	 else closedir(dir_);
	 std::string casename_rmap = output_dir + filenameWithoutExtension;

    std::string header_file = casename_dem + ".header";
    std::ifstream dem_input(INPUT_DEM);
    std::ofstream header_output(header_file);

    if (!dem_input.is_open() || !header_output.is_open()) {
        std::cerr << "Error creating header file." << std::endl;
        MPI_Finalize();		  
        return 1;
    }

    for (int i = 0; i < 6; i++) {
        std::getline(dem_input, line);
        header_output << line << std::endl;
    }

    dem_input.close();
    header_output.close();

    int nlines = nrows / NFILES;
    int rem = nrows % NFILES;
    std::vector<int> line_numbers(NFILES);
    int sum = 0;

    for (int i = 0; i < NFILES - 1; i++) {
        line_numbers[i] = sum + nlines;
        sum += nlines;

        if (i < rem && rem > 0) {
            line_numbers[i]++;
            sum++;
        }
    }
	 line_numbers[NFILES - 1] = nrows;

	int num_lines = line_numbers.size();
	int lines_per_process = (num_lines + size - 1) / size;
	int start_idx = rank * lines_per_process;
	int end_idx = std::min((rank + 1) * lines_per_process, num_lines);
	
	if(rank==size-1){ //last file
		end_idx=num_lines;
	}

	split_dem_to_bin(casename_dem, start_idx, end_idx, line_numbers, ncols );
	if (IS_MANN == "YES") {
   	split_mann_to_bin(casename_mann, start_idx, end_idx, line_numbers, ncols);
	}
   if (IS_RMAP == "YES") {
   	split_rmap_to_bin(casename_rmap, start_idx, end_idx, line_numbers, ncols);
	}

	

	 MPI_Barrier(MPI_COMM_WORLD);


	 if(rank==0){
    	std::cout << "ASCII files generated" << std::endl;
	 }
	 
	 MPI_Finalize();


    return 0;
}

