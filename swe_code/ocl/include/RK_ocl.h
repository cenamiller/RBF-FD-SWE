#ifndef RK_OCL_H
#define RK_OCL_H

#include <profiling.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <swe_config.h>
#include <device_setup.h>

// Use kernels to calculate the Runge-Kutta K value and add together to solve the differential equation
void RK_substep_ocl(cl_kernel, cl_kernel, cl_kernel, cl_kernel, cl_kernel, cl_command_queue, int);

// Calculate righthand side of the differential equation
void eval_RHS_ocl(cl_kernel, cl_command_queue);

void copy_arr_ocl(cl_kernel, cl_command_queue);

// Adds each RK K to D- D stores weighted sum of K values for each timestep
void update_D_ocl(cl_kernel, cl_command_queue, const fType);

// Calculate the RK K value
void eval_K_ocl(cl_kernel, cl_command_queue, const fType);

// Add weighted sum of K values to H
void update_H_ocl(cl_kernel, cl_command_queue);

#endif
