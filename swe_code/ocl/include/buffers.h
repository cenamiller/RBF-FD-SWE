#ifndef BUFFERS_H
#define BUFFERS_H

#include <device_setup.h>
#include <stdlib.h>
#include <string.h>
#include <swe_config.h>

// Creates buffer buff with contents of buff_source
void load_buffer_fType(cl_context context, cl_command_queue commandQueue, cl_mem* buff, fType* buff_source, int buff_size);

void load_buffer_int(cl_context context, cl_command_queue commandQueue, cl_mem* buff, int* buff_source, int buff_size);

void load_all_buffers(cl_context context, cl_command_queue commandQueue, LPSMD_buffers* LPSMD_buffs, cl_mem* F_buff, cl_mem* D_buff);

void release_all_buffers(LPSMD_buffers* LPSMD_buffs, cl_mem* H_buff, cl_mem* F_buff, cl_mem* D_buff, cl_mem* K_buff);
#endif
