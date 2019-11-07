// Written by Samuel Elliott, Summer 2017.

#include <rcm.h>

#include <stdlib.h>
#include <stdio.h>

// returns integer id array pointer to the rcm ordered mapping
int* rcm_mapping(GSMD_struct GSMD) {

    // root node id for BFS
    int r = 0;

    // relevant atm struct constants
    int Nnodes = GSMD.Nnodes;
    int Nnbr = GSMD.Nnbr;

    // DM stencil mapping
    int* idx = GSMD.idx;

    int* mapping = (int *) malloc(sizeof(int)*Nnodes);
    int* mapped = (int *) calloc(Nnodes,sizeof(int));
    int* degree = (int *) malloc(sizeof(int)*Nnbr);
    int* x = (int *) malloc(sizeof(int)*Nnbr);
	
    // step through until all nodes are mapped
    mapping[0] = r;
    mapped[r] = 1;
    int mapping_iter = 1;

    for (int i = 0; i < Nnodes; i++) {
        int node_id = mapping[i];
        
        find_degree(idx, mapped, degree, x, node_id, Nnbr);
        reorder_increasing(x, degree, Nnbr, Nnodes);

        for (int inbr = 1; inbr < Nnbr; inbr++) {
            int nbr_id = x[inbr];
            if (mapped[nbr_id] == 0) {
                mapping[mapping_iter] = nbr_id;
                mapped[nbr_id] = 1;
                mapping_iter++;
                if (mapping_iter == Nnodes) break;
            }
        }

        if (mapping_iter == Nnodes) break;
    }

    // reverse ordering
    for (int i = 0; i < (Nnodes-1)/2; i++) {
        int temp_id = mapping[i];
        mapping[i] = mapping[Nnodes-1-i];
        mapping[Nnodes-1-i] = temp_id;
    }

    // print_mapping(mapping,Nnodes);
    free(x);
    free(degree);
    free(mapped);

    return mapping;
}


// prints info on the maximum bandwidth of the RCM reordered sparse DM matrix
// note that this function should be called with a remapped idx
void print_max_bandwidth(GSMD_struct GSMD) {

    int Nnodes = GSMD.Nnodes;
    int Nnbr = GSMD.Nnbr;

    // DM stencil mapping
    int* idx = GSMD.idx;

    int max_bw = 0;
    int sum = 0;

    for (int i = 0; i < Nnodes; i++) {
        int max_id = 0, min_id = Nnodes;
        for (int j = 0; j < Nnbr; j++) {
            int id = idx[(i * Nnbr) + j];
            if (id > max_id) max_id = id;
            if (id < min_id) min_id = id;
        }

        int bw = max_id - min_id;
        sum += bw;
        if (bw > max_bw) max_bw = bw;
    }

    printf("\n\n=========================== RCM MATRIX INFO ======================\
            \n\nMAX BANDWIDTH = %d\nAVERAGE_BANDWIDTH=%d\
            \n\n==================================================================\n\n"
            ,max_bw,sum/Nnodes);
}

void print_mapping(int* mapping, int Nnodes) {

    printf("\n");
    for (int i = 0; i < Nnodes; i++) {
        printf("%6d  --------> %6d\n", i, mapping[i]);
    }
    printf("\n\n");
}

void find_degree(const int* idx, const int* mapped, int* degree, int* x, const int node_id, const int Nnbr) {

    const int* y1 = &idx[Nnbr * node_id];
    for (int inbr = 0; inbr < Nnbr; inbr++) {
        degree[inbr] = 0;
        x[inbr] = y1[inbr];
    }

    for (int inbr = 1; inbr < Nnbr; inbr++) {
        int nbr_node_id = x[inbr];
        const int* y2 = &idx[Nnbr * nbr_node_id];
        for (int jnbr = 1; jnbr < Nnbr; jnbr++) if (mapped[y2[jnbr]] == 0) degree[inbr]++;
    }
}

void min_max(int* degree, int* min, int* max, int Nnbr, int offset) {

    min[0] = Nnbr + 1;
    max[0] = 0;

    for (int i = offset; i < Nnbr; i++) {
        if (degree[i] < min[0]) min[0] = degree[i];
        if (degree[i] > max[0]) max[0] = degree[i];
    }
}

void reorder_increasing(int* x, int* degree, const int Nnbr, const int Nnodes) {

    int iter = 1;
    int min, max, dummy_id;
    for (int i = 1; i < Nnbr; i++) {
        min_max(degree, &min, &max, Nnbr, i);
        for (int j = i; j < Nnbr; j++) {
            if (degree[j] == min) {
                dummy_id = x[iter];
                x[iter] = x[j];
                x[j] = dummy_id;
                degree[j] = degree[iter];
                degree[iter] = min;
                iter++;
                break;
            }
        }
    }
	
    int count;
    for (int i = 1; i < Nnbr; i+=count) {
        min_max(degree, &min, &max, Nnbr, i);
        count = 0;
        for (int j = 1; j < Nnbr; j++) {
            if (degree[j] == min) count++;
        }

        for (int j = i; j < i+count; j++) {
            int min_id = Nnodes, min_k = Nnbr;
            for (int k = j; k < i+count; k++) {
                if (x[k] < min_id) {
                    min_id = x[k];
                    min_k = k;
                }
            }
            dummy_id = x[j];
            x[j] = x[min_k];
            x[min_k] = dummy_id;
        }
    }
}
