// Written by Richelle Streater, Summer 2018- openCL version of rk4_rbffd_swe.c in main folder.

#ifdef USE_OCL
#include <RK_ocl.h>

extern PSMD_struct* LPSMD;
extern timing_struct local_timer;

void RK_substep_ocl(cl_kernel eval_RHS_kernel, cl_kernel copy_arr_kernel, cl_kernel update_D_kernel, cl_kernel eval_K_kernel, cl_kernel update_H_kernel, cl_command_queue commandQueue, int substep_id) {
    
    // full RK4 timestep length
    const fType dt = LPSMD->dt;
    
    switch(substep_id) {
        case 0:
            // get F_0 = d/dt(K_0 = H)
            eval_RHS_ocl(eval_RHS_kernel, commandQueue);
            // initialize D = F_0
            copy_arr_ocl(copy_arr_kernel, commandQueue);
            // evaluate K_1 = H + (dt/2) * F_0
            eval_K_ocl(eval_K_kernel, commandQueue, dt/2.0);
            break;
            
        case 1:
            // get F_1 = d/dt(K_1 = H + (dt/2) * F_0)
            eval_RHS_ocl(eval_RHS_kernel, commandQueue);
            // update D += 2 * F_1
            update_D_ocl(update_D_kernel, commandQueue, 2.0);
            // evaluate K_2 = H + (dt/2) * F_0
            eval_K_ocl(eval_K_kernel, commandQueue, dt/2.0);
            break;
            
        case 2:
            // get F_2 = d/dt(K_2 = H + (dt/2) * F_1)
            eval_RHS_ocl(eval_RHS_kernel, commandQueue);
            // update D += 2 * F_2
            update_D_ocl(update_D_kernel, commandQueue, 2.0);
            // evaluate K_3 = H + dt * F_0
            eval_K_ocl(eval_K_kernel, commandQueue, dt);
            break;
            
        case 3:
            // get F_3 = d/dt(K_3 = H + dt * F_2)
            eval_RHS_ocl(eval_RHS_kernel, commandQueue);
            // update D += F_3
            update_D_ocl(update_D_kernel, commandQueue, 1.0);
            // add everything together to get H_n+1
            update_H_ocl(update_H_kernel, commandQueue);
            break;
    }
}

void eval_RHS_ocl(cl_kernel kernel, cl_command_queue commandQueue){
    cl_int err = 0;
    double t_start = getTime();
    
    // The CFDL=1, SFDL=0 kernel is tiled differently than the other kernels
    #if ((SIMD_LENGTH != 4) && (SIMD_LENGTH != 8)) || (defined(USE_CFDL) && !defined(USE_SFDL))
    int tile_length = 1;
    
    #else
    int tile_length = SIMD_LENGTH;
    #endif
    
    // global size and work-group size- wgsize is flexible within device limits.
    size_t wgsize = 128;
    size_t globalSize = (((LPSMD->compute_size)/tile_length)/wgsize + 1) * wgsize;
    
    // Launch kernel
    err = clEnqueueNDRangeKernel(commandQueue, kernel, 1, NULL, &globalSize, &wgsize, 0, NULL, NULL);
    
    if (err != CL_SUCCESS){
        printf("Error in clEnqueueNDRangeKernel: %s\n", checkError(err));
        exit(-1);
    }
    commandQueueSync(commandQueue);
    local_timer.t_eval_rhs[local_timer.attempt] +=  (getTime() - t_start);
}

void copy_arr_ocl(cl_kernel kernel, cl_command_queue commandQueue){
    cl_int err = 0;
    
    // global size and work-group size
    size_t wgsize = 64;
    size_t globalSize = (((LPSMD->compute_size)*4)/ wgsize + 1) * wgsize;
    
    // Launch kernel
    err = clEnqueueNDRangeKernel(commandQueue, kernel, 1, NULL, &globalSize, &wgsize, 0, NULL, NULL);
    
    if (err != CL_SUCCESS){
        printf("Error in clEnqueueNDRangeKernel: %s\n", checkError(err));
        exit(-1);
    }
    commandQueueSync(commandQueue);
}

void update_D_ocl(cl_kernel kernel, cl_command_queue commandQueue, const fType coefficient){
    cl_int err = 0;
    double t_start = getTime();
    
    // global size and work-group size
    size_t wgsize = 64;
    size_t globalSize = (((LPSMD->compute_size)*4)/ wgsize + 1) * wgsize;
    
    err = clSetKernelArg(kernel, 2, sizeof(fType), (void*)&coefficient);
    
    // Launch kernel
    err = clEnqueueNDRangeKernel(commandQueue, kernel, 1, NULL, &globalSize, &wgsize, 0, NULL, NULL);
    
    if (err != CL_SUCCESS){
        printf("Error in clEnqueueNDRangeKernel: %s\n", checkError(err));
        exit(-1);
    }
    commandQueueSync(commandQueue);
    local_timer.t_update_D[local_timer.attempt] += (getTime() - t_start);
}

void eval_K_ocl(cl_kernel kernel, cl_command_queue commandQueue, const fType dt){
    cl_int err = 0;
    double t_start = getTime();
    
    // global size and work-group size
    size_t wgsize = 64;
    size_t globalSize = (((LPSMD->compute_size)*4)/ wgsize + 1) * wgsize;
    
    err = clSetKernelArg(kernel, 3, sizeof(fType), (void*)&dt);
    
    // Launch kernel
    err = clEnqueueNDRangeKernel(commandQueue, kernel, 1, NULL, &globalSize, &wgsize, 0, NULL, NULL);
  
    if (err != CL_SUCCESS){
        printf("Error in clEnqueueNDRangeKernel: %s\n", checkError(err));
        exit(-1);
    }
    commandQueueSync(commandQueue);
    local_timer.t_eval_K[local_timer.attempt] +=  (getTime() - t_start);
}

void update_H_ocl(cl_kernel kernel, cl_command_queue commandQueue){
    cl_int err = 0;
    double t_start = getTime();
    
    // global size and work-group size
    size_t wgsize = 64;
    size_t globalSize = (((LPSMD->compute_size)*4)/wgsize + 1) * wgsize;
    
    err = clEnqueueNDRangeKernel(commandQueue, kernel, 1, NULL, &globalSize, &wgsize, 0, NULL, NULL);
    if (err != CL_SUCCESS){
        printf("Error in clEnqueueNDRangeKernel: %s\n", checkError(err));
        exit(-1);
    }
    
    commandQueueSync(commandQueue);
    local_timer.t_update_H[local_timer.attempt] +=  (getTime() - t_start);
}
#endif
