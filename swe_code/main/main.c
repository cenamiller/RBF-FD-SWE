// Written by Samuel Elliott, Summer 2017. Last update- Richelle Streater, 8/29/18

#include <swe_config.h>
#include <runtime_params.h>
#include <input.h>
#include <rcm.h>
#include <reorder_nodes.h>
#include <layout.h>
#include <matrix_transformations.h>
#include <init_patches.h>
#include <rk4_rbffd_swe.h>
#include <halos.h>
#include <profiling.h>

#ifdef USE_OCL
#include <device_setup.h>
#include <buffers.h>
#include <RK_ocl.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef _OPENACC
#include <openacc.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Global Variables
sim_params_struct sim_params;		// runtime parameterizations and model configuration options
PSMD_struct* LPSMD;		        	// local patch static model data
timing_struct local_timer;	    	// timing data for proviling


// mpi rank and comm size
int mpi_rank;
int mpi_size;

void verify_output(fType* H1, fType* H2);

int main(int argc, char** argv) {
	
    // ------------------------------ MPI initialization ------------------------------------
    
    #ifdef USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    #else
    mpi_rank = 0;
    mpi_size = 1;
    #endif
    
    double t_start = getTime();
    double t_start_inner;
	
    // --------------------------------------------------------------------------------------

    // ------------------ Get Simulation and Model Parameterizations ------------------------
    // initialize simulation parameterization struct
    get_runtime_params();

    // initialize local mpi timer
    local_timer = init_timer(sim_params.nattempts);
	
    if (mpi_rank == 0) print_sim_params();

    #ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    // ----------------------------------------------------------------------------------- //

    // =================================================================================== //
    //                                                                                     //
    //                          Initialize Simulation Datasets                             //
    //                                                                                     //
    // =================================================================================== //

    fType* H_init;
    fType* H_100;
    int* mapping;
    fType* H_init_global;
    fType* H_100_global;
    GSMD_struct GSMD;
    LPSMD = (PSMD_struct*) malloc(sizeof(PSMD_struct));
    
    // =========================== Rank 0 Initial Operations ============================= //

    if (mpi_rank == 0) {

        // ------------------------ Read Simulation Input ------------------------------------ //

        // read atmospheric data from the input file
        if (sim_params.USE_NETCDF == 1) {
            
            #ifdef USE_NCIO
            GSMD = get_GSMD_nc((char*) &sim_params.inputFile[0]);
            H_init_global = get_ICs_nc((char*) &sim_params.inputFile[0]);
            H_100_global = get_FCs_nc((char*) &sim_params.inputFile[0]);
            #else
            printf("Executable was not built for NetCDF I/O compatibility.\n\n");
            exit(0);
            #endif
        }
        else {
            GSMD = get_GSMD_bin((char*) &sim_params.inputFile[0]);
            H_init_global = get_ICs_bin((char*) &sim_params.inputFile[0]);
            H_100_global = get_FCs_bin((char*) &sim_params.inputFile[0]);
        }

        // --------------------------------------------------------------------------------------

        // ------------------------ Perform RCM Node Reordering ---------------------------------

        // array of node ids of the corresponding rcm node mapping
        if (sim_params.USE_RCM == 1)
        {
            // get rcm ordering of nodeset
            mapping = rcm_mapping(GSMD);

            // reorder all necessary variables for the rcm node ordering
            GSMD = reorder_GSMD_struct(mapping, GSMD);
            H_init_global = reorder_2D_fp_arr(H_init_global, mapping, GSMD.Nnodes, 4);
            H_100_global = reorder_2D_fp_arr(H_100_global, mapping, GSMD.Nnodes, 4);

        }
        
        // --------------------------------------------------------------------------------------

        // ----------------- Perform DM and SVM Layout Transformations --------------------------

        // perform any necessary padding
        GSMD = pad_GSMD_data(GSMD);
        H_init_global = pad_2D_fp_arr(H_init_global, GSMD.Nnodes, GSMD.padded_Nnodes, 4, 4);
        H_100_global = pad_2D_fp_arr(H_100_global, GSMD.Nnodes, GSMD.padded_Nnodes, 4, 4);

        //If using SFDL then perform transposition/tiling of all necessary data
        #ifdef USE_SFDL
        GSMD = convert2_SFDL(GSMD);
        #endif

        // If not using CFDL then transpose/tile H_init_global
        // note that tiling is expected to have negligible performance effects but makes MPI simpler
        #ifndef USE_CFDL
        H_init_global = transpose_2D_fp_arr(H_init_global, GSMD.padded_Nnodes, 4);
        H_init_global = tile_2D_fp_arr(H_init_global, 4, GSMD.padded_Nnodes, SIMD_LENGTH);
        H_100_global = transpose_2D_fp_arr(H_100_global, GSMD.padded_Nnodes, 4);
        H_100_global = tile_2D_fp_arr(H_100_global, 4, GSMD.padded_Nnodes, SIMD_LENGTH);
        #endif
        
        // --------------------------------------------------------------------------------------

        // ---------------------------- MPI Partitioning ----------------------------------------
		
        // get global patch decomposition (GPD)
        patch_struct* GPD = get_GPD(GSMD);
        verify_valid_GPD(GPD);
        print_patches(GPD);

        // initialize local PSMD
        init_LPSMD(GPD, GSMD);

        // initialize H_init
        H_init = get_patch_ICs(H_init_global, GPD);
        H_100 = get_patch_ICs(H_100_global, GPD);

        // --------------------------------------------------------------------------------------

    }
    else {
		
        // ----------------------- Initialize MPI Patch Data ------------------------------------

        // initialize local PSMD
        patch_struct* dummy_GPD;
        init_LPSMD(dummy_GPD, GSMD);

        H_init = get_patch_ICs(H_init_global, dummy_GPD);
        H_100 = get_patch_ICs(H_100_global, dummy_GPD);
		
    }

    // --------------------------------------------------------------------------------------

    local_timer.t_init = getTime() - t_start;

    // -------------------------------- Timestepping ----------------------------------------
    
    // initialize timestepping vars
    fType* H = (fType*) malloc(sizeof(fType) * LPSMD->patch_size * 4);
    fType* K = (fType*) malloc(sizeof(fType) * LPSMD->patch_size * 4);
    fType* F = (fType*) malloc(sizeof(fType) * LPSMD->compute_size * 4);
    fType* D = (fType*) malloc(sizeof(fType) * LPSMD->compute_size * 4);

    #ifdef _OPENACC
		
    // Copy in all data from LPSMD struct to device
    #pragma acc enter data copyin(LPSMD[:1])
    #pragma acc enter data copyin(\
        LPSMD->x[:LPSMD->compute_size],LPSMD->y[:LPSMD->compute_size],LPSMD->z[:LPSMD->compute_size],LPSMD->f[:LPSMD->compute_size],LPSMD->ghm[:LPSMD->compute_size],\
        LPSMD->gradghm[:LPSMD->compute_size*3],LPSMD->p_u[:LPSMD->compute_size*3],LPSMD->p_v[:LPSMD->compute_size*3],LPSMD->p_w[:LPSMD->compute_size*3],\
        LPSMD->idx[:LPSMD->compute_size*LPSMD->padded_Nnbr],LPSMD->Dx[:LPSMD->compute_size*LPSMD->padded_Nnbr],LPSMD->Dy[:LPSMD->compute_size*LPSMD->padded_Nnbr],\
        LPSMD->Dz[:LPSMD->compute_size*LPSMD->padded_Nnbr],LPSMD->L[:LPSMD->compute_size*LPSMD->padded_Nnbr])

    // State variable and related data
    #pragma acc enter data copyin(H[:LPSMD->patch_size*4],K[:LPSMD->patch_size*4],F[:LPSMD->compute_size*4],D[:LPSMD->compute_size*4],H_init[:LPSMD->patch_size*4])
    #endif
    
    // -------------------------- Setup OpenCL Kernels/Environment ----------------------------

    #ifdef USE_OCL
    
    // OpenCL environment
    cl_context context;
    cl_command_queue commandQueue;
    cl_platform_id platform;
    cl_uint deviceCount;
    cl_device_id device;
    cl_program program;
    
    // Define kernels
    cl_kernel eval_RHS_kernel;
    cl_kernel copy_arr_kernel;
    cl_kernel update_D_kernel;
    cl_kernel eval_K_kernel;
    cl_kernel update_H_kernel;
    
    // Free all opencl buffers
    LPSMD_buffers* LPSMD_buffs = (LPSMD_buffers*) malloc(sizeof(LPSMD_buffers));
    cl_mem H_buff;
    cl_mem F_buff;
    cl_mem K_buff;
    cl_mem D_buff;

    fType* H_cl;
    fType* K_cl;
    
    // Use the first device and platform that the program finds
    int deviceID = 0;
    int platformID = 0;
    
    #endif
    
    if (sim_params.OCL == 1) {
        
        #ifdef USE_OCL
       
        // Time the OpenCL initialization
        t_start = getTime();
    
        createPlatform(&platform, platformID);
        createContext(&context, platform, &device, deviceID);
        createCommandQueue(context, device, &commandQueue);
    
        // Find kernels
        char kernelFile[MAX_FILE_PATH_SIZE];
        if (getenv("SWE_KERNEL_FILE") == NULL) {
            strcpy(kernelFile, getenv("SWE_TOP_DIR"));
            #if (defined(USE_CFDL) && !defined(USE_SFDL))
            strcat(kernelFile, "/swe_code/ocl/kernel_CFDL.cl");
            #else
            strcat(kernelFile, "/swe_code/ocl/kernel.cl");
            #endif
        }
        else strcpy(kernelFile,getenv("SWE_KERNEL_FILE"));

        // Compile all kernels
        compileKernel(kernelFile, "eval_RHS_kernel", context, &program, &device, &eval_RHS_kernel);
        compileKernel(kernelFile, "copy_arr_kernel", context, &program, &device, &copy_arr_kernel);
        compileKernel(kernelFile, "update_D_kernel", context, &program, &device, &update_D_kernel);
        compileKernel(kernelFile, "eval_K_kernel", context, &program, &device, &eval_K_kernel);
        compileKernel(kernelFile, "update_H_kernel", context, &program, &device, &update_H_kernel);
    
        load_all_buffers(context, commandQueue, LPSMD_buffs, &F_buff, &D_buff);
        local_timer.t_ocl = getTime() - t_start;

        #else
        // Check before running that the program was built with OpenCL
        if (sim_params.OCL == 1) {
            printf("Executable was not built for OpenCL compatibility.\n\n");
            exit(0);
        }
        #endif
        
    }
    
    // --------------------------------------------------------------------------------------
    
    for (int attempt = 0; attempt < sim_params.nattempts; attempt++) {
        
        // reset initial state varible data for each attempt
        copy_fp_arr(H, H_init, LPSMD->patch_size * 4);
      
        if (sim_params.OCL == 1) {
            
            #ifdef USE_OCL
      
            // Reset H and K vectors- leave maps open for MPI implementation
            H_buff = createBuffer(context, CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR, sizeof(fType) * LPSMD->patch_size*4);
            H_cl = (fType*)mapBuffer(commandQueue, H_buff, CL_MAP_WRITE, sizeof(fType) * LPSMD->patch_size*4);
            memcpy(H_cl, H, sizeof(fType) * LPSMD->patch_size*4);

            K_buff = createBuffer(context, CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR, sizeof(fType) * LPSMD->patch_size*4);
            K_cl = (fType*)mapBuffer(commandQueue, K_buff, CL_MAP_WRITE, sizeof(fType) * LPSMD->patch_size*4);
            memcpy(K_cl, H, sizeof(fType) * LPSMD->patch_size*4);
            
            if (mpi_size <= 1) {
                unmapBuffer(commandQueue, H_buff, H_cl);
                unmapBuffer(commandQueue, K_buff, K_cl);
            }
            
            // Set kernel arguments
            set_eval_RHS_args(eval_RHS_kernel, LPSMD_buffs, K_buff, F_buff, LPSMD->gh0, LPSMD->compute_size, LPSMD->padded_Nnbr, LPSMD->Nnbr, LPSMD->compute_pid_s);
            set_copy_arr_args(copy_arr_kernel, D_buff, F_buff, 4*LPSMD->compute_size);
            set_update_D_args(update_D_kernel, F_buff, D_buff, LPSMD->compute_size);
            set_eval_K_args(eval_K_kernel, H_buff, F_buff, LPSMD->compute_size, K_buff, LPSMD->compute_pid_s);
            set_update_H_args(update_H_kernel, H_buff, D_buff, K_buff, LPSMD->compute_size, LPSMD->dt, LPSMD->compute_pid_s);
            #endif
            
        }
        
        #ifdef USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
        #endif
        
        local_timer.attempt = attempt;
        t_start = getTime();

        // ------------------------------ Runge Kutta Main Loop ---------------------------------
        
        for (int i = 1; i <= sim_params.nsteps; i++) {
            for (int rk_val = 0; rk_val < 4; rk_val++) {
                if (sim_params.OCL == 1) {
                    
                    #ifdef USE_OCL
                    RK_substep_ocl(eval_RHS_kernel, copy_arr_kernel, update_D_kernel, eval_K_kernel, update_H_kernel, commandQueue, rk_val);
                    
                    #ifdef MPI
                    if (mpi_size > 1) {
                    
                        // Temporary MPI implementation with OCL- Use mapping to exchange SVM halo's
                        t_start_inner = getTime();
                        exchange_SVM_halos(K_cl);
                        if (rk_val == 3) exchange_SVM_halos(H_cl);
                        local_timer.t_mpi[local_timer.attempt] += getTime() - t_start_inner;
                    }
                    #endif
                    #endif
                    
                }
                else RK_substep(H, K, F, D, rk_val);
            }
        }
        
        // --------------------------------------------------------------------------------------
        
        if (sim_params.OCL == 1) { 
            
            #ifdef USE_OCL
            if (mpi_size > 1){
                unmapBuffer(commandQueue, H_buff, H_cl);
                unmapBuffer(commandQueue, K_buff, K_cl);
            }
            #endif
        }
        
        local_timer.t_main[local_timer.attempt] = getTime() - t_start;
        
        if (sim_params.OCL == 1) {
            
            #ifdef USE_OCL
            
            // Copy contents of H_buff back to H array
            H = (fType*)mapBuffer(commandQueue, H_buff, CL_MAP_READ, sizeof(fType)*LPSMD->patch_size*4);
            unmapBuffer(commandQueue, H_buff, H);
            #endif
            
        }

        #ifdef _OPENACC
        #pragma acc update self(H[:LPSMD->patch_size*4])
        #endif
        
        // Compare H array to expected values at t=100
        if (mpi_rank == 0) printf("Attempt %d:\t", attempt + 1);
        verify_output(H, H_100);
    }
    
    if (sim_params.OCL == 1) {

        #ifdef USE_OCL
   
        // Release OpenCL objects
        if(commandQueue)clReleaseCommandQueue(commandQueue);
        if(context)clReleaseContext(context);
        if(program)clReleaseProgram(program);
    
        if(eval_RHS_kernel)clReleaseKernel(eval_RHS_kernel);
        if(copy_arr_kernel)clReleaseKernel(copy_arr_kernel);
        if(update_D_kernel)clReleaseKernel(update_D_kernel);
        if(eval_K_kernel)clReleaseKernel(eval_K_kernel);
        if(update_H_kernel)clReleaseKernel(update_H_kernel);
    
        release_all_buffers(LPSMD_buffs, &H_buff, &F_buff, &D_buff, &K_buff);
        free(LPSMD_buffs);
        #endif
    }

    process_profiling_results();

    // --------------------------------------------------------------------------------------

    // ------------------------ Cleanup and MPI finalization --------------------------------
	
    #ifdef _OPENACC
    #pragma acc exit data delete(H[:LPSMD->patch_size*4],K[:LPSMD->patch_size*4],F[:LPSMD->compute_size*4],D[:LPSMD->compute_size*4],H_init[:LPSMD->patch_size*4])
    #pragma acc exit data delete(\
        LPSMD->x[:LPSMD->compute_size],LPSMD->y[:LPSMD->compute_size],LPSMD->z[:LPSMD->compute_size],LPSMD->f[:LPSMD->compute_size],LPSMD->ghm[:LPSMD->compute_size],\
        LPSMD->gradghm[:LPSMD->compute_size*3],LPSMD->p_u[:LPSMD->compute_size*3],LPSMD->p_v[:LPSMD->compute_size*3],LPSMD->p_w[:LPSMD->compute_size*3],\
        LPSMD->idx[:LPSMD->compute_size*LPSMD->padded_Nnbr],LPSMD->Dx[:LPSMD->compute_size*LPSMD->padded_Nnbr],LPSMD->Dy[:LPSMD->compute_size*LPSMD->padded_Nnbr],\
        LPSMD->Dz[:LPSMD->compute_size*LPSMD->padded_Nnbr],LPSMD->L[:LPSMD->compute_size*LPSMD->padded_Nnbr])
    #endif

    #ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    #endif

    // --------------------------------------------------------------------------------------

}

// Check output H with known value at t=100
void verify_output(fType* H1, fType* H2) {

    int count = 0;
    int SV_IDS[4] = {U_SV_ID, V_SV_ID, W_SV_ID, H_SV_ID};
    for (int i = LPSMD->compute_pid_s; i < LPSMD->compute_pid_e; i++) {
        for (int j = 0; j < 4; j++) {

            // get absolute value of error
            fType diff = H1[SVM_LIN_ID(i,SV_IDS[j])] - H2[SVM_LIN_ID(i,SV_IDS[j])];
            diff = diff < 0 ? -diff : diff;
			
            // ensure error is acceptable (and check for nan!)
            if (diff > TOLERANCE || diff != diff ) {
                count++;
            }
        }
    }

    int global_count;

    #ifdef USE_MPI
    MPI_Reduce((const void *) &count, (void*) &global_count, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    #else
    global_count = count;
    #endif

    if (mpi_rank == 0) {
        if (global_count != 0) {
        printf("\n\nVERIFICATION FAILED: There were a total of %d values that exceeded an error of %.1e\n\n", global_count, TOLERANCE);
        }
        else printf("VERIFICATION SUCCEEDED\n");
    }

    fflush(stdout);
    
    #ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
}
