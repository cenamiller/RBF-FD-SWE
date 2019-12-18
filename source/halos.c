// Written by Samuel Elliott, Summer 2017.

#include <halos.h>

#ifdef USE_MPI
#include <mpi.h>

// local patch static model data
extern PSMD_struct* LPSMD;

void exchange_SVM_halos(fType* H) {

    MPI_Request request[4];
    MPI_Status status[4];

    #ifdef _OPENACC
    #pragma acc host_data use_device(H)
    #endif
    {
    
    // send neighbor halo layer data
    // left neighbors halo layer
    MPI_Isend((void*) &H[LPSMD->compute_pid_s*4], LPSMD->lnbr_halo_size*4, MPI_FTYPE, LPSMD->lnbr_rank, 0, MPI_COMM_WORLD, &request[0]);
    // right neighbors halo layer
    MPI_Isend((void*) &H[LPSMD->rnbr_halo_pid_s*4], LPSMD->rnbr_halo_size*4, MPI_FTYPE, LPSMD->rnbr_rank, 1, MPI_COMM_WORLD, &request[1]);

    // initialize halo layer MPI recieves
    // right halo layer
    MPI_Irecv((void*) &H[LPSMD->compute_pid_e*4], LPSMD->rhalo_size*4, MPI_FTYPE, LPSMD->rnbr_rank, 0, MPI_COMM_WORLD, &request[2]);
    // left halo layer
    MPI_Irecv((void*) &H[LPSMD->halo_pid_s*4], LPSMD->lhalo_size*4, MPI_FTYPE, LPSMD->lnbr_rank, 1, MPI_COMM_WORLD, &request[3]);
    }

    MPI_Waitall(4, request, status);
}

#endif
