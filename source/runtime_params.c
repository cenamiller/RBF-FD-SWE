// Written by Samuel Elliott, Summer 2017

#include <runtime_params.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

extern sim_params_struct sim_params;	// parameterizations and model configuration options

// gets runtime simulation/model parameterizations from runtime environment
void get_runtime_params() {
	
    // input file
    strcpy(&sim_params.inputFile[0],getenv("SWE_INPUT_FILE"));
	
    // IO method
    sim_params.USE_NETCDF = atoi(getenv("SWE_USE_NETCDF"));
    
    // How many devices to split into
    sim_params.num_compute_units = atoi(getenv("SWE_NUM_COMPUTE_UNITS"));
	
    // output file
    #ifdef USE_HIST
    strcpy(&sim_params.outputFile[0],getenv("SWE_OUTPUT_FILE"));
    sim_params.hist_int = atoi(getenv("SWE_HIST_INT"));
    #endif

    // number of simulation steps
    sim_params.nsteps = atoi(getenv("SWE_NSTEPS"));

    // number of simulation attempts
    sim_params.nattempts = atoi(getenv("SWE_NATTEMPTS"));

    // node ordering and layout options
    sim_params.USE_RCM = atoi(getenv("SWE_USE_RCM"));

    // Use opencl for parallelization
    sim_params.OCL = atoi(getenv("SWE_USE_OCL"));
}

// print runtime parameterizations at beggining of simulation
void print_sim_params() {

    printf("\n\n=============================================== Runtime Configurations ===================================================\n\n");
    printf("Input from file: yes\n\tInput File Path: %s\n\n",&sim_params.inputFile[0]);

    // IO method
    #ifdef USE_NCIO
    printf("I/O Method: %s\n",sim_params.USE_NETCDF == 1 ? "NETCDF" : "BINARY");
    #endif

    // output file info
    #ifdef USE_HIST
    printf("Write history to output file: yes\n\tOutput File Path: %s\n\tHistory Interval (timesteps): %d\n\n",&sim_params.outputFile[0],sim_params.hist_int);
    #endif

    // number of simulation steps
    printf("Simulation length (timesteps): %d\n",sim_params.nsteps);
    printf("Number of attempts: %d\n",sim_params.nattempts);

    // node ordering and layout options
    printf("USE RCM Node Ordering: %s\n", sim_params.USE_RCM == 1 ? "yes" : "no");
    printf("\n============================================================================================================================\n\n");
}
