#ifndef DEVICE_SETUP_H
#define DEVICE_SETUP_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <swe_config.h>

void createPlatform(cl_platform_id*, int);

void createContext(cl_context*, cl_platform_id, cl_device_id*, int);

void createCommandQueue(cl_context, cl_device_id, cl_command_queue*);

void compileKernel(const char*, const char*, cl_context, cl_program*, cl_device_id*, cl_kernel*);

cl_mem createBuffer(cl_context, cl_mem_flags, size_t);

void* mapBuffer(cl_command_queue, cl_mem, cl_map_flags, size_t);

void unmapBuffer(cl_command_queue, cl_mem, void*);

void commandQueueSync(cl_command_queue);

// Read contents of file as cstring, and add preprocessor directives to kernel cstring
char* readKernelSource(const char*);

const char* checkError(cl_int);

void set_eval_RHS_args(cl_kernel, LPSMD_buffers*, cl_mem, cl_mem, fType, int, int, int, int);

void set_copy_arr_args(cl_kernel, cl_mem, cl_mem, const int);

void set_update_D_args(cl_kernel, cl_mem, cl_mem, const int);

void set_eval_K_args(cl_kernel, cl_mem, cl_mem, const int, cl_mem, int);

void set_update_H_args (cl_kernel, cl_mem, cl_mem, cl_mem, const int, const fType, int);

#endif
