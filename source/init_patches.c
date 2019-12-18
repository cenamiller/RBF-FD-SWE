// Written by Samuel Elliott, Summer 2017. Last updated- Richelle Streater, Summer 2018.

#include <init_patches.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

// mpi rank and comm size
extern int mpi_size; //(1 in serial case)
extern int mpi_rank; //(0 in serial case)

// local patch static model data
extern PSMD_struct* LPSMD;

// determines/assigns global patch decomposition (GPD) for MPI
patch_struct* get_GPD(GSMD_struct GSMD) {

    // allocate array of patch_structs for each mpi rank
    patch_struct* GPD = (patch_struct*) malloc(sizeof(patch_struct) * mpi_size);

    // ======================================== Determine Compute Patch Distributions ============================================== //

    int g_size = GSMD.padded_Nnodes;							                // global domain size
    int p_size = ROUNDDOWN_SIMD(g_size / mpi_size);						        // default patch size
    int n_default = mpi_size - ((g_size - (p_size * mpi_size)) / SIMD_LENGTH);  // number of patches with default size

    // set compute gids for patches with default size
    for (int rank = 0; rank < n_default; rank++) {
        GPD[rank].compute_gid_s = p_size * rank;
        GPD[rank].compute_gid_e = p_size * (rank + 1);
        GPD[rank].compute_size = p_size;
    }
    // set compute gids for remaining patches (default size + SIMD_LENGTH)
    for (int rank = n_default; rank < mpi_size; rank++) {
        GPD[rank].compute_gid_s = (p_size * n_default) + ((p_size + SIMD_LENGTH) * (rank - n_default));
        GPD[rank].compute_gid_e = (p_size * n_default) + ((p_size + SIMD_LENGTH) * ((rank + 1) - n_default));
        GPD[rank].compute_size = p_size + SIMD_LENGTH;
    }

    // ============================================================================================================================= //

    // ================================================ Determine Halo Layers ====================================================== //
	
    // for each patches compute domain find the smallest and largest gids needed by the DMs to determine the halo layer extents
    for (int rank = 0; rank < mpi_size; rank++) {
	
        // initialize halo layer gids to that of the compute domain
        int min_gid = GPD[rank].compute_gid_s;
        int max_gid = GPD[rank].compute_gid_e - 1;

        // iterate through each compute node
        for (int gid = GPD[rank].compute_gid_s; gid < GPD[rank].compute_gid_e; gid ++) {

            // iterate through each neighbor
            for (int nbr = 0; nbr < GSMD.Nnbr; nbr++) {

                // get gid of neighbor
                int nbr_gid = GSMD.idx[DM_LIN_ID(gid, nbr, GSMD.padded_Nnbr)];

                // if new max/min id update values
                if (nbr_gid > max_gid) max_gid = nbr_gid;
                if (nbr_gid < min_gid) min_gid = nbr_gid;
            }
        }
        // use resulting min and max gids for the patch halo layer
        // note we extend the halo layer to SIMD boundaries
        GPD[rank].halo_gid_s = ROUNDDOWN_SIMD(min_gid);
        GPD[rank].halo_gid_e = ROUNDUP_SIMD(max_gid + 1);

        // total patch size
        GPD[rank].patch_size = GPD[rank].halo_gid_e - GPD[rank].halo_gid_s;
    }
    
    // update neighbor halo layer info
    for (int rank = 0; rank < mpi_size; rank++) {

        // left neighbor halo layer info
        if (rank == 0) GPD[rank].lnbr_halo_gid_e = 0;
        else GPD[rank].lnbr_halo_gid_e = GPD[rank - 1].halo_gid_e;

        // right neighbor halo layer info
        if (rank == mpi_size - 1) GPD[rank].rnbr_halo_gid_s = GSMD.padded_Nnodes;
		else GPD[rank].rnbr_halo_gid_s = GPD[rank + 1].halo_gid_s;
    }

    return GPD;
}

void init_LPSMD(patch_struct* GPD, GSMD_struct GSMD) {
    #ifdef USE_MPI
    
    // MPI nonblocking comm vars
    MPI_Request* request;
    MPI_Status* status;
    int n_send_recv;

    // ========================================= Get Global Constants and Paraneters ==================================== //

    n_send_recv = mpi_rank == 0 ? mpi_size - 1 : 1;

    request = (MPI_Request*) malloc(sizeof(MPI_Request) * n_send_recv);
    status = (MPI_Status*) malloc(sizeof(MPI_Status) * n_send_recv);
    #endif
	
    if (mpi_rank == 0) {

        // assign rank 0 constants/parameters
        LPSMD->gh0 = GSMD.gh0;
        LPSMD->dt = GSMD.dt;
        LPSMD->Nnbr = GSMD.Nnbr;
        LPSMD->padded_Nnbr = GSMD.padded_Nnbr;

        #ifdef USE_MPI
        
        // send constants/parameters to all ranks
        for (int rank = 1; rank < mpi_size; rank++) {
            MPI_Isend((void*) &LPSMD[0], sizeof(PSMD_struct), MPI_BYTE, rank, rank, MPI_COMM_WORLD,&request[rank - 1]);
        }
        #endif
    }
    else {
        #ifdef USE_MPI
        MPI_Irecv((void*) &LPSMD[0], sizeof(PSMD_struct), MPI_BYTE, 0, mpi_rank, MPI_COMM_WORLD, request);
        #endif
    }
    
    #ifdef USE_MPI
    MPI_Waitall(n_send_recv, request, status);
    free(request);
    free(status);
    #endif
    
    // ====================================== Get MPI Patch Decomposition Info ============================================== //

    // local mpi patch decomposition info within global context
    patch_struct local_patch;
    #ifdef USE_MPI
    n_send_recv = mpi_rank == 0 ? mpi_size + 1 : 1;

    request = (MPI_Request*) malloc(sizeof(MPI_Request) * n_send_recv);
    status = (MPI_Status*) malloc(sizeof(MPI_Status) * n_send_recv);

    // send/recieve patch docomposition info
    if (mpi_rank == 0) {

        // send patch information to all other ranks
        for (int rank = 0; rank < mpi_size; rank++) {
            MPI_Isend((void*) &GPD[rank], sizeof(patch_struct), MPI_BYTE, rank, rank, MPI_COMM_WORLD,&request[rank + 1]);
        }
    }
	
    MPI_Irecv((void*) &local_patch, sizeof(patch_struct), MPI_BYTE, 0, mpi_rank, MPI_COMM_WORLD, request);
    MPI_Waitall(n_send_recv, request, status);
    free(request);
    free(status);
	
    // ======================================== Set Local PSMD MPI Data =================================================== //

    // MPI neighbor rank data
    LPSMD->lnbr_rank = mpi_rank - 1;
    LPSMD->rnbr_rank = mpi_rank + 1;

    // first/last ranks have no left/right neighbors
    if (mpi_rank == 0) LPSMD->lnbr_rank = MPI_PROC_NULL;
	if (mpi_rank == mpi_size - 1) LPSMD->rnbr_rank = MPI_PROC_NULL;

    #else
    LPSMD->lnbr_rank = mpi_rank; // RDL: there are no neighbors in the serial case.
    LPSMD->rnbr_rank = mpi_rank;
    local_patch = GPD[0]; // RDL: I think this might emulate the message passing that is going on.
    #endif

    // use local_patch to determine patch info
    // note that we choose compute_pid_s = 0 for local patch indexing

    // sizes
    LPSMD->patch_size = local_patch.patch_size;
    LPSMD->compute_size = local_patch.compute_gid_e - local_patch.compute_gid_s;
    LPSMD->lhalo_size = local_patch.compute_gid_s - local_patch.halo_gid_s;
    LPSMD->rhalo_size = local_patch.halo_gid_e - local_patch.compute_gid_e;
    LPSMD->lnbr_halo_size = local_patch.lnbr_halo_gid_e - local_patch.compute_gid_s;
    LPSMD->rnbr_halo_size = local_patch.compute_gid_e - local_patch.rnbr_halo_gid_s;

    // convert gids to pids using compute_pid_s = 0 for local patch indexing
    LPSMD->halo_pid_s = 0;
    LPSMD->compute_pid_s = LPSMD->halo_pid_s + LPSMD->lhalo_size;
    LPSMD->compute_pid_e = LPSMD->compute_pid_s + LPSMD->compute_size;
    LPSMD->halo_pid_e = LPSMD->compute_pid_e + LPSMD->rhalo_size;
    LPSMD->rnbr_halo_pid_s = LPSMD->compute_pid_e - LPSMD->rnbr_halo_size;
    LPSMD->lnbr_halo_pid_e = LPSMD->compute_pid_s + LPSMD->lnbr_halo_size;

    // ======================================== Allocate space for Local PSMD Data ====================================== //

    // ---------------------------- static vars ---------------------------------- //

    // scalars
    LPSMD->x = (fType*) malloc(sizeof(fType) * LPSMD->compute_size);
    LPSMD->y = (fType*) malloc(sizeof(fType) * LPSMD->compute_size);
    LPSMD->z = (fType*) malloc(sizeof(fType) * LPSMD->compute_size);
    LPSMD->f = (fType*) malloc(sizeof(fType) * LPSMD->compute_size);
    LPSMD->ghm = (fType*) malloc(sizeof(fType) * LPSMD->compute_size);

    // 3D vars
    LPSMD->gradghm = (fType*) malloc(3 * sizeof(fType) * LPSMD->compute_size);
    LPSMD->p_u = (fType*) malloc(3 * sizeof(fType) * LPSMD->compute_size);
    LPSMD->p_v = (fType*) malloc(3 * sizeof(fType) * LPSMD->compute_size);
    LPSMD->p_w = (fType*) malloc(3 * sizeof(fType) * LPSMD->compute_size);

    // ---------------------------- DM data -------------------------------------- //

    LPSMD->idx = (int*) malloc(sizeof(int) * LPSMD->compute_size * LPSMD->padded_Nnbr);
    LPSMD->Dx = (fType*) malloc(sizeof(fType) * LPSMD->compute_size * LPSMD->padded_Nnbr);
    LPSMD->Dy = (fType*) malloc(sizeof(fType) * LPSMD->compute_size * LPSMD->padded_Nnbr);
    LPSMD->Dz = (fType*) malloc(sizeof(fType) * LPSMD->compute_size * LPSMD->padded_Nnbr);
    LPSMD->L = (fType*) malloc(sizeof(fType) * LPSMD->compute_size * LPSMD->padded_Nnbr);

    // ========================================== Send/Recv Local PSMD Data ============================================= //

    #ifdef USE_MPI
    int n_datasets = 14;
    n_send_recv = mpi_rank == 0 ? (mpi_size + 1) * n_datasets : n_datasets;

    request = (MPI_Request*) malloc(sizeof(MPI_Request) * n_send_recv);
    status = (MPI_Status*) malloc(sizeof(MPI_Status) * n_send_recv);

    // send/recieve patch docomposition info
    if (mpi_rank == 0) {
        
        // send patch information to all other ranks
        for (int rank = 0; rank < mpi_size; rank++) {
            patch_struct patch = GPD[rank];
            MPI_Isend((void*) &GSMD.x[patch.compute_gid_s * 1], patch.compute_size * 1, MPI_FTYPE, rank, (mpi_size * 0) + rank, MPI_COMM_WORLD, &request[((mpi_size + 1) * 0) + rank + 1]);
            MPI_Isend((void*) &GSMD.y[patch.compute_gid_s * 1], patch.compute_size * 1, MPI_FTYPE, rank, (mpi_size * 1) + rank, MPI_COMM_WORLD, &request[((mpi_size + 1) * 1) + rank + 1]);
            MPI_Isend((void*) &GSMD.z[patch.compute_gid_s * 1], patch.compute_size * 1, MPI_FTYPE, rank, (mpi_size * 2) + rank, MPI_COMM_WORLD, &request[((mpi_size + 1) * 2) + rank + 1]);
            MPI_Isend((void*) &GSMD.f[patch.compute_gid_s * 1], patch.compute_size * 1, MPI_FTYPE, rank, (mpi_size * 3) + rank, MPI_COMM_WORLD, &request[((mpi_size + 1) * 3) + rank + 1]);
            MPI_Isend((void*) &GSMD.ghm[patch.compute_gid_s * 1], patch.compute_size * 1, MPI_FTYPE, rank, (mpi_size * 4) + rank, MPI_COMM_WORLD, &request[((mpi_size + 1) * 4) + rank + 1]);
			
            MPI_Isend((void*) &GSMD.p_u[patch.compute_gid_s * 3], patch.compute_size * 3, MPI_FTYPE, rank, (mpi_size * 5) + rank, MPI_COMM_WORLD, &request[((mpi_size + 1) * 5) + rank + 1]);
            MPI_Isend((void*) &GSMD.p_v[patch.compute_gid_s * 3], patch.compute_size * 3, MPI_FTYPE, rank, (mpi_size * 6) + rank, MPI_COMM_WORLD, &request[((mpi_size + 1) * 6) + rank + 1]);
            MPI_Isend((void*) &GSMD.p_w[patch.compute_gid_s * 3], patch.compute_size * 3, MPI_FTYPE, rank, (mpi_size * 7) + rank, MPI_COMM_WORLD, &request[((mpi_size + 1) * 7) + rank + 1]);
            MPI_Isend((void*) &GSMD.gradghm[patch.compute_gid_s * 3], patch.compute_size * 3, MPI_FTYPE, rank, (mpi_size * 8) + rank, MPI_COMM_WORLD, &request[((mpi_size + 1) * 8) + rank + 1]);

            MPI_Isend((void*) &GSMD.idx[patch.compute_gid_s * GSMD.padded_Nnbr], patch.compute_size * GSMD.padded_Nnbr, MPI_INT, rank, (mpi_size * 9) + rank, MPI_COMM_WORLD, &request[((mpi_size + 1) * 9) + rank + 1]);
            MPI_Isend((void*) &GSMD.Dx[patch.compute_gid_s * GSMD.padded_Nnbr], patch.compute_size * GSMD.padded_Nnbr, MPI_FTYPE, rank, (mpi_size * 10) + rank, MPI_COMM_WORLD, &request[((mpi_size + 1) * 10) + rank + 1]);
            MPI_Isend((void*) &GSMD.Dy[patch.compute_gid_s * GSMD.padded_Nnbr], patch.compute_size * GSMD.padded_Nnbr, MPI_FTYPE, rank, (mpi_size * 11) + rank, MPI_COMM_WORLD, &request[((mpi_size + 1) * 11) + rank + 1]);
            MPI_Isend((void*) &GSMD.Dz[patch.compute_gid_s * GSMD.padded_Nnbr], patch.compute_size * GSMD.padded_Nnbr, MPI_FTYPE, rank, (mpi_size * 12) + rank, MPI_COMM_WORLD, &request[((mpi_size + 1) * 12) + rank + 1]);
            MPI_Isend((void*) &GSMD.L[patch.compute_gid_s * GSMD.padded_Nnbr], patch.compute_size * GSMD.padded_Nnbr, MPI_FTYPE, rank, (mpi_size * 13) + rank, MPI_COMM_WORLD, &request[((mpi_size + 1) * 13) + rank + 1]);
        }
    }

    MPI_Irecv((void*) LPSMD->x, LPSMD->compute_size * 1, MPI_FTYPE, 0, (mpi_size * 0) + mpi_rank, MPI_COMM_WORLD, &request[(n_send_recv / n_datasets) * 0]);
    MPI_Irecv((void*) LPSMD->y, LPSMD->compute_size * 1, MPI_FTYPE, 0, (mpi_size * 1) + mpi_rank, MPI_COMM_WORLD, &request[(n_send_recv / n_datasets) * 1]);
    MPI_Irecv((void*) LPSMD->z, LPSMD->compute_size * 1, MPI_FTYPE, 0, (mpi_size * 2) + mpi_rank, MPI_COMM_WORLD, &request[(n_send_recv / n_datasets) * 2]);
    MPI_Irecv((void*) LPSMD->f, LPSMD->compute_size * 1, MPI_FTYPE, 0, (mpi_size * 3) + mpi_rank, MPI_COMM_WORLD, &request[(n_send_recv / n_datasets) * 3]);
    MPI_Irecv((void*) LPSMD->ghm, LPSMD->compute_size * 1, MPI_FTYPE, 0, (mpi_size * 4) + mpi_rank, MPI_COMM_WORLD, &request[(n_send_recv / n_datasets) * 4]);

    MPI_Irecv((void*) LPSMD->p_u, LPSMD->compute_size * 3, MPI_FTYPE, 0, (mpi_size * 5) + mpi_rank, MPI_COMM_WORLD, &request[(n_send_recv / n_datasets) * 5]);
    MPI_Irecv((void*) LPSMD->p_v, LPSMD->compute_size * 3, MPI_FTYPE, 0, (mpi_size * 6) + mpi_rank, MPI_COMM_WORLD, &request[(n_send_recv / n_datasets) * 6]);
    MPI_Irecv((void*) LPSMD->p_w, LPSMD->compute_size * 3, MPI_FTYPE, 0, (mpi_size * 7) + mpi_rank, MPI_COMM_WORLD, &request[(n_send_recv / n_datasets) * 7]);
    MPI_Irecv((void*) LPSMD->gradghm, LPSMD->compute_size * 3, MPI_FTYPE, 0, (mpi_size * 8) + mpi_rank, MPI_COMM_WORLD, &request[(n_send_recv / n_datasets) * 8]);
    
    MPI_Irecv((void*) LPSMD->idx, LPSMD->compute_size * LPSMD->padded_Nnbr, MPI_INT, 0, (mpi_size * 9) + mpi_rank, MPI_COMM_WORLD, &request[(n_send_recv / n_datasets) * 9]);
    MPI_Irecv((void*) LPSMD->Dx, LPSMD->compute_size * LPSMD->padded_Nnbr, MPI_FTYPE, 0, (mpi_size * 10) + mpi_rank, MPI_COMM_WORLD, &request[(n_send_recv / n_datasets) * 10]);
    MPI_Irecv((void*) LPSMD->Dy, LPSMD->compute_size * LPSMD->padded_Nnbr, MPI_FTYPE, 0, (mpi_size * 11) + mpi_rank, MPI_COMM_WORLD, &request[(n_send_recv / n_datasets) * 11]);
    MPI_Irecv((void*) LPSMD->Dz, LPSMD->compute_size * LPSMD->padded_Nnbr, MPI_FTYPE, 0, (mpi_size * 12) + mpi_rank, MPI_COMM_WORLD, &request[(n_send_recv / n_datasets) * 12]);
    MPI_Irecv((void*) LPSMD->L, LPSMD->compute_size * LPSMD->padded_Nnbr, MPI_FTYPE, 0, (mpi_size * 13) + mpi_rank, MPI_COMM_WORLD, &request[(n_send_recv / n_datasets) * 13]);

    MPI_Waitall(n_send_recv, request, status);
    free(request);
    free(status);
    
    // ==================================================================== Modify idx for local patch ids ================================================================= //

    for (int i = 0; i < LPSMD->compute_size; i++) {
        for (int j = 0; j < LPSMD->padded_Nnbr; j++) {
            LPSMD->idx[(i * LPSMD->padded_Nnbr) + j] = LPSMD->idx[(i * LPSMD->padded_Nnbr) + j] - local_patch.halo_gid_s;
        }
    }
    
    #else
    
    // Copy all data members from GSMD to LPSMD in non-MPI case
    memcpy(LPSMD->x, GSMD.x, sizeof(fType) * LPSMD->compute_size);
    memcpy(LPSMD->y, GSMD.y, sizeof(fType) * LPSMD->compute_size);
    memcpy(LPSMD->z, GSMD.z, sizeof(fType) * LPSMD->compute_size);
    memcpy(LPSMD->f, GSMD.f, sizeof(fType) * LPSMD->compute_size);
    memcpy(LPSMD->ghm, GSMD.ghm, sizeof(fType) * LPSMD->compute_size);

    memcpy(LPSMD->p_u, GSMD.p_u, 3 * sizeof(fType) * LPSMD->compute_size);
    memcpy(LPSMD->p_v, GSMD.p_v, 3 * sizeof(fType) * LPSMD->compute_size);
    memcpy(LPSMD->p_w, GSMD.p_w, 3 * sizeof(fType) * LPSMD->compute_size);
    memcpy(LPSMD->gradghm, GSMD.gradghm, 3 * sizeof(fType) * LPSMD->compute_size);
	
    memcpy(LPSMD->idx, GSMD.idx, sizeof(int) * LPSMD->compute_size * LPSMD->padded_Nnbr);
    memcpy(LPSMD->Dx, GSMD.Dx, sizeof(fType) * LPSMD->compute_size * LPSMD->padded_Nnbr);
    memcpy(LPSMD->Dy, GSMD.Dy, sizeof(fType) * LPSMD->compute_size * LPSMD->padded_Nnbr);
    memcpy(LPSMD->Dz, GSMD.Dz, sizeof(fType) * LPSMD->compute_size * LPSMD->padded_Nnbr);
    memcpy(LPSMD->L, GSMD.L, sizeof(fType) * LPSMD->compute_size * LPSMD->padded_Nnbr);

    #endif
}

fType* get_patch_ICs(fType* H_global, patch_struct* GPD) {

    // allocate array for dataset
    fType* H_local = (fType*) malloc(sizeof(fType) * LPSMD->patch_size * 4);

    #ifdef USE_MPI
    
    // MPI nonblocking comm vars
    MPI_Request* request;
    MPI_Status* status;
    int n_send_recv;
    
    n_send_recv = mpi_rank == 0 ? mpi_size + 1 : 1;

    request = (MPI_Request*) malloc(sizeof(MPI_Request) * n_send_recv);
    status = (MPI_Status*) malloc(sizeof(MPI_Status) * n_send_recv);

    if (mpi_rank == 0) {

        // send patch information to all other ranks
        for (int rank = 0; rank < mpi_size; rank++) {
            patch_struct patch = GPD[rank];
            MPI_Isend((void*) &H_global[patch.halo_gid_s * 4], patch.patch_size * 4, MPI_FTYPE, rank, rank, MPI_COMM_WORLD, &request[rank + 1]);
        }
    }

    MPI_Irecv((void*) H_local, LPSMD->patch_size * 4, MPI_FTYPE, 0, mpi_rank, MPI_COMM_WORLD, &request[0]);

    MPI_Waitall(n_send_recv, request, status);
    free(request);
    free(status);

    #else
    memcpy(H_local, H_global, sizeof(fType) * LPSMD->patch_size * 4);
    #endif
	
    return H_local;
}

void verify_valid_GPD(patch_struct* GPD) {

    for (int rank = 0; rank < mpi_size; rank++) {
        if (GPD[rank].lnbr_halo_gid_e > GPD[rank].compute_gid_e || GPD[rank].rnbr_halo_gid_s < GPD[rank].compute_gid_s) {
            printf("\n\nOverdecoposition Error: Too many MPI processes. Halo layers overlap between patches.\n\n");
            fflush(stdout);
            
            #ifdef USE_MPI
            MPI_Abort(MPI_COMM_WORLD, 1);
            #endif
        }
    }
}

void print_patches(patch_struct* GPD) {

    printf("\n\n====================================================================================== MPI Patch Decomposition ================================================================================================\n");

    for (int rank = 0; rank < mpi_size; rank++) {
        printf("\n\tRank #%3d:\thalo_gid_s = %7d,\tcompute_gid_s = %7d,\tlnbr_halo_gid_e = %7d,\trnbr_halo_gid_s = %7d,\tcompute_gid_e = %7d,\t\thalo_gid_e = %7d", rank, GPD[rank].halo_gid_s, GPD[rank].compute_gid_s, GPD[rank].lnbr_halo_gid_e, GPD[rank].rnbr_halo_gid_s, GPD[rank].compute_gid_e, GPD[rank].halo_gid_e);
    }

    printf("\n\n================================================================================================================================================================================================================\n\n");
    return;
}
