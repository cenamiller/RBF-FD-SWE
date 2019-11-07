// Written by Richelle Streater, Summer 2018.

#ifdef USE_OCL
#include <buffers.h>

extern PSMD_struct* LPSMD;

// Copy contents of array into buffer by creating a new array, mapping the buffer to the new array, and copying source array
// Read array of fTypes
void load_buffer_fType(cl_context context, cl_command_queue commandQueue, cl_mem* buff, fType* buff_source, int buff_size) {
    *buff = createBuffer(context, CL_MEM_READ_ONLY|CL_MEM_ALLOC_HOST_PTR, sizeof(fType) * buff_size);
    fType* buff_cl = (fType*)mapBuffer(commandQueue, *buff, CL_MAP_WRITE, sizeof(fType) * buff_size);
    memcpy(buff_cl, buff_source, sizeof(fType) * buff_size);
    unmapBuffer(commandQueue, *buff, buff_cl);
}

// Read array of ints
void load_buffer_int(cl_context context, cl_command_queue commandQueue, cl_mem *buff, int* buff_source, int buff_size) {
    *buff = createBuffer(context, CL_MEM_READ_ONLY|CL_MEM_ALLOC_HOST_PTR, sizeof(int) * buff_size);
    int* buff_cl = (int*)mapBuffer(commandQueue, *buff, CL_MAP_WRITE, sizeof(int) * buff_size);
    memcpy(buff_cl, buff_source, sizeof(int) * buff_size);
    unmapBuffer(commandQueue, *buff, buff_cl);
}

// Copy input arrays to buffers and initialize the RK kernel inputs
void load_all_buffers(cl_context context, cl_command_queue commandQueue, LPSMD_buffers* LPSMD_buffs, cl_mem* F_buff, cl_mem* D_buff) {
    int padded_Nnodes = LPSMD->compute_size;
    int padded_Nnbr = LPSMD->padded_Nnbr;
    int Nnbr = LPSMD->Nnbr;
    
    load_buffer_int(context, commandQueue, &LPSMD_buffs->idx, LPSMD->idx, padded_Nnodes * padded_Nnbr);
    load_buffer_fType(context, commandQueue, &LPSMD_buffs->Dx, LPSMD->Dx, padded_Nnodes * padded_Nnbr);
    load_buffer_fType(context, commandQueue, &LPSMD_buffs->Dy, LPSMD->Dy, padded_Nnodes * padded_Nnbr);
    load_buffer_fType(context, commandQueue, &LPSMD_buffs->Dz, LPSMD->Dz, padded_Nnodes * padded_Nnbr);
    load_buffer_fType(context, commandQueue, &LPSMD_buffs->L, LPSMD->L, padded_Nnodes * padded_Nnbr);
    load_buffer_fType(context, commandQueue, &LPSMD_buffs->x, LPSMD->x, padded_Nnodes);
    load_buffer_fType(context, commandQueue, &LPSMD_buffs->y, LPSMD->y, padded_Nnodes);
    load_buffer_fType(context, commandQueue, &LPSMD_buffs->z, LPSMD->z, padded_Nnodes);
    load_buffer_fType(context, commandQueue, &LPSMD_buffs->f, LPSMD->f, padded_Nnodes);
    load_buffer_fType(context, commandQueue, &LPSMD_buffs->ghm, LPSMD->ghm, padded_Nnodes);
    load_buffer_fType(context, commandQueue, &LPSMD_buffs->p_u, LPSMD->p_u, padded_Nnodes * 3);
    load_buffer_fType(context, commandQueue, &LPSMD_buffs->p_v, LPSMD->p_v, padded_Nnodes * 3);
    load_buffer_fType(context, commandQueue, &LPSMD_buffs->p_w, LPSMD->p_w, padded_Nnodes * 3);
    load_buffer_fType(context, commandQueue, &LPSMD_buffs->gradghm, LPSMD->gradghm, padded_Nnodes * 3);
    
    // F and D start out as NULL
    *F_buff = createBuffer(context, CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR, sizeof(fType)*padded_Nnodes*4);
    *D_buff = createBuffer(context, CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR, sizeof(fType)*padded_Nnodes*4);
}

// Free memory for all buffers
void release_all_buffers(LPSMD_buffers* LPSMD_buffs, cl_mem* H_buff, cl_mem* F_buff, cl_mem* D_buff, cl_mem* K_buff) {
    
    // For each buffer, check if it exists, and if not, release
    if(LPSMD_buffs->x)clReleaseMemObject(LPSMD_buffs->x);
    if(LPSMD_buffs->y)clReleaseMemObject(LPSMD_buffs->y);
    if(LPSMD_buffs->z)clReleaseMemObject(LPSMD_buffs->z);
    if(LPSMD_buffs->f)clReleaseMemObject(LPSMD_buffs->f);
    if(LPSMD_buffs->ghm)clReleaseMemObject(LPSMD_buffs->ghm);
    if(LPSMD_buffs->p_u)clReleaseMemObject(LPSMD_buffs->p_u);
    if(LPSMD_buffs->p_v)clReleaseMemObject(LPSMD_buffs->p_v);
    if(LPSMD_buffs->p_w)clReleaseMemObject(LPSMD_buffs->p_w);
    if(LPSMD_buffs->gradghm)clReleaseMemObject(LPSMD_buffs->gradghm);
    if(LPSMD_buffs->idx)clReleaseMemObject(LPSMD_buffs->idx);
    if(LPSMD_buffs->Dx)clReleaseMemObject(LPSMD_buffs->Dx);
    if(LPSMD_buffs->Dy)clReleaseMemObject(LPSMD_buffs->Dy);
    if(LPSMD_buffs->Dz)clReleaseMemObject(LPSMD_buffs->Dz);
    if(LPSMD_buffs->L)clReleaseMemObject(LPSMD_buffs->L);
    
    if(*F_buff)clReleaseMemObject(*F_buff);
    if(*K_buff)clReleaseMemObject(*K_buff);
    if(*D_buff)clReleaseMemObject(*D_buff);
    if(*H_buff)clReleaseMemObject(*H_buff);
}

#endif
