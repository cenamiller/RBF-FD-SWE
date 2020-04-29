

// Authorship:
//
// Written by Samuel Elliott, Summer 2017.
//
// Refactored for more generality by Rich Loft, November, 2019

#include <stdlib.h>
#include <stdio.h>

#include <rcm.h>

// returns integer id array pointer to the rcm ordered mapping

int* rcm_mapping(const int NNodes, const int NNbrs_Padded, const int NNbrs, int *idx) {

  // root node id for BFS
  int r = 0;
  int* mapping = (int *) malloc(sizeof(int)*NNodes);
  int* mapped = (int *) calloc(NNodes,sizeof(int));
  int* degree = (int *) malloc(sizeof(int)*NNbrs);
  int* x = (int *) malloc(sizeof(int)*NNbrs);
  
  // step through until all nodes are mapped
  mapping[0] = r;
  mapped[r] = 1;
  int mapping_iter = 1;

  for (int i = 0; i < NNodes; i++) {
    int node_id = mapping[i];
        
    rcm_find_degree(idx, mapped, degree, x, node_id, NNbrs_Padded, NNbrs);
    rcm_reorder_increasing(x, degree, NNbrs, NNodes);

    for (int inbr = 1; inbr < NNbrs; inbr++) {
      int nbr_id = x[inbr];
      if (mapped[nbr_id] == 0) {
	mapping[mapping_iter] = nbr_id;
	mapped[nbr_id] = 1;
	mapping_iter++;
	if (mapping_iter == NNodes) break;
      }
    }

    if (mapping_iter == NNodes) break;
  }

  // reverse ordering
  for (int i = 0; i < (NNodes-1)/2; i++) {
    int temp_id = mapping[i];
    mapping[i] = mapping[NNodes-1-i];
    mapping[NNodes-1-i] = temp_id;
  }

  // print_mapping(mapping,NNodes);
  free(x);
  free(degree);
  free(mapped);

  return mapping;
}


// prints info on the maximum bandwidth of the RCM reordered sparse DM matrix
// note that this function should be called with a remapped idx
void rcm_print_max_bandwidth(const int NNodes, const int NNbrs_Padded, const int NNbrs, int *idx) {

  int max_bw = 0;
  int sum = 0;

  for (int i = 0; i < NNodes; i++) {
    int max_id = 0, min_id = NNodes;
    for (int j = 0; j < NNbrs; j++) {
      int id = idx[(i * NNbrs_Padded) + j];
      if (id > max_id) max_id = id;
      if (id < min_id) min_id = id;
    }

    int bw = max_id - min_id;
    sum += bw;
    if (bw > max_bw) max_bw = bw;
  }

  printf("\nMATRIX INFO: MAX BANDWIDTH = %d AVERAGE_BANDWIDTH=%d\n\n"
	 ,max_bw,sum/NNodes);
}

void rcm_print_mapping(int* mapping, int NNodes) {

  printf("\n");
  for (int i = 0; i < NNodes; i++) {
    printf("%6d  --------> %6d\n", i, mapping[i]);
  }
  printf("\n\n");
}

// Check RCM mapping for errors

void rcm_check_mapping(int* mapping, int NNodes) {

  int *one_to_one =   (int *) calloc(NNodes, sizeof(int));

  for (int i = 0; i < NNodes; i++) {
    if (mapping[i] >= NNodes || mapping[i] < 0){
      printf("RCM ERROR: Node remapping [%i] out of bounds. Exiting...\n",i);
      exit(-1);
    }
    else{
      one_to_one[mapping[i]] += 1;
    }
  }

  for (int i = 0; i < NNodes; i++) {
    if (one_to_one[i] != 1){
      printf("RCM ERROR: Node remapping [%i] not one to one. Exiting...\n",i);
      exit(-1);
    }
  }

  printf("RCM mapping check passed...\n");

  free(one_to_one);

}

void rcm_find_degree(const int* idx, const int* mapped, int* degree, int* x, const int node_id, const int NNbrs_Padded, const int NNbrs) {

  const int* y1 = &idx[NNbrs_Padded * node_id];
  for (int inbr = 0; inbr < NNbrs; inbr++) {
    degree[inbr] = 0;
    x[inbr] = y1[inbr];
  }

  for (int inbr = 1; inbr < NNbrs; inbr++) {
    int nbr_node_id = x[inbr];
    const int* y2 = &idx[NNbrs_Padded * nbr_node_id];
    for (int jnbr = 1; jnbr < NNbrs; jnbr++){
      if (mapped[y2[jnbr]] == 0){
	degree[inbr]++;
      }
    }

  }
}

void rcm_min_max(int* degree, int* min, int* max, int NNbrs, int offset) {

  min[0] = NNbrs + 1;
  max[0] = 0;

  for (int i = offset; i < NNbrs; i++) {
    if (degree[i] < min[0]) min[0] = degree[i];
    if (degree[i] > max[0]) max[0] = degree[i];
  }
}

void rcm_reorder_increasing(int* x, int* degree, const int NNbrs, const int NNodes) {

  int iter = 1;
  int min, max, dummy_id;

  for (int i = 1; i < NNbrs; i++) {
    rcm_min_max(degree, &min, &max, NNbrs, i);
    for (int j = i; j < NNbrs; j++) {
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
  for (int i = 1; i < NNbrs; i+=count) {
    rcm_min_max(degree, &min, &max, NNbrs, i);
    count = 0;
    for (int j = 1; j < NNbrs; j++) {
      if (degree[j] == min) count++;
    }

    for (int j = i; j < i+count; j++) {
      int min_id = NNodes, min_k = NNbrs;
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
