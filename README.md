**Two-dimensional Runoff Inundation Toolkit for Operational Needs (TRITON; Morales-Hernández et al., Submitted to Environmental Modelling & Software)**

[TRITON Website](https://triton.ornl.gov/)

*Refer to User's Guide, located in doc subdirectory, for instructions on requirements, installation, and other details.* 

A 2D open source flood simulation tool designed for modern high performance computing (HPC). The core of the tool is a computationally efficient, physics-based hydraulic model that operates on a regular/structured grid and solves the full 2D shallow water equations. The key features of TRITON are:


*  It can operate on multiple computer platforms and utilize modern HPC environments. The users can take advantage of:
1.  Implementation with a single central processing unit (CPU) or multiple CPUs (using OpenMP+MPI)
2.  Implementation with a single graphics processing unit (GPU) or multiple GPUs (using CUDA+MPI)

Highest TRITON computational efficiency can be achieved by using GPU implementation.


*  TRITON utilizes topographical data (e.g., digital elevation model [DEM], light detection and ranging [LIDAR]), as its base input, in a uniform (Cartesian) grid structure. The model can be driven by streamflow hydrographs at specified locations or gridded runoff hydrographs, or both which serves as the model’s hydrological forcing. The primary TRITON output includes water depth and 2D velocity maps at user defined time intervals. Other variables such as unit discharge values can be outputted. TRITON can also output timeseries of simulated results as user-defined point locations.

*  TRITON is developed on Linux/Unix platform. The input/output files can be either in ascii or binary formats. A set of tools and instruction are provided for format conversion.

*  The model utilizes International System of Units (SI). Users who are more familiar with United States (US) customary units need to perform proper unit conversion themselves.


*  A set of Test Cases are provided :

1.  Case 01 : TaumSauk
2.  Case 02a: Paraboloid Case - Resolution 0.04m
3.  Case 02b: Paraboloid Case - Resolution 0.02m
4.  Case 02c: Paraboloid Case - Resolution 0.01m
5.  Case 02d: Paraboloid Case - Resolution 0.005m
6.  Case 03 : Runoff
7.  Case 04	: Harvey 30m 
8.  Case 05	: Harvey 10m


Morales-Hernández, M., Sharif, M.B., Kalyanapu, A., Ghafoor, S.K., Dullo, T.T., Gangrade, S., Kao, S.C., Norman, M.R. and Evans, K.J., 2021. TRITON: A Multi-GPU Open Source 2D Hydrodynamic Flood Model. Environmental Modelling & Software, p.105034 [https://doi.org/10.1016/j.envsoft.2021.105034](https://doi.org/10.1016/j.envsoft.2021.105034)


**Instructions to run a case with parallel input**

1.  Open the script "scriptSplitASCII" and configure the following parameters (in caps):

*  TRITON_DIR: triton directory
*  NFILES: The number of files to split the ASCII file(s) into. It should match the number of ranks in the TRITON simulation
*  INPUT_DEM: DEM input file (in ASCII format)
*  IS_MANN: flag for the mann file (YES or NO)
*  INPUT_MANN: in case IS_MANN=YES, this parameter contains the mann input file (in ASCII format)
*  IS_RMAP: flag for the rmap file (YES or NO)
*  INPUT_RMAP: in case IS_RMAP=YES, this parameter contains the rmap input file (in ASCII format)
*  OUTPUT_FORMAT: ASC or BIN, depending on the desired output format for the sequence of files
*  ASCII2BIN_FOLDER: in case OUTPUT_FORMAT=BIN, it points to the ascii2bin directory (not the ascii2bin_dem)
*  ASCII2BIN_RMAP_FOLDER: in case OUTPUT_FORMAT=BIN and IS_RMAP=YES, it points to the ascii2bin_rmap directory



2.  Run the script.  


3.  Configure the cfg file as follows for the following parameters:

*  dem_filename: include the path of the split files plus the basename (without underscores or numbers). Example: "input/dem/bin/par/case03"
*  Add header_filename and point it to the exact name of the header that the script has generated. Example: "input/dem/bin/par/case03.header" 
*  n_infile: include the path of the split files plus the basename (without underscores or numbers). Example: "input/mann/bin/par/case03"
*  runoff_map: include the path of the split files plus the basename (without underscores or numbers). Example: ""input/runoff/bin/par/case03_runoff"
*  input_format: BIN
*  Add input_option: PAR

4.  Run TRITON with the same MPI ranks as stated in NFILES.


