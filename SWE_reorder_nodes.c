#include <stdlib.h>
#include <SWE.h>
#include <reorder_nodes.h>

void SWE_reorder_nodes(SWE_struct *SWE, int *mapping){

  // replace all variable data pointed to by the GSMD struct with reordered versions
  // scalar variables

  SWE->x = reorder_1D_fp_arr(SWE->x, mapping, SWE->NNodes);
  SWE->y = reorder_1D_fp_arr(SWE->y, mapping, SWE->NNodes);
  SWE->z = reorder_1D_fp_arr(SWE->z, mapping, SWE->NNodes);
  SWE->f = reorder_1D_fp_arr(SWE->f, mapping, SWE->NNodes);
  SWE->ghm = reorder_1D_fp_arr(SWE->ghm, mapping, SWE->NNodes);

  // vector variables
  SWE->p_u     = reorder_2D_fp_arr(SWE->p_u, mapping, SWE->NNodes, 3, 3);
  SWE->p_v     = reorder_2D_fp_arr(SWE->p_v, mapping, SWE->NNodes, 3, 3);
  SWE->p_w     = reorder_2D_fp_arr(SWE->p_w, mapping, SWE->NNodes, 3, 3);
  SWE->gradghm = reorder_2D_fp_arr(SWE->gradghm, mapping, SWE->NNodes, 3, 3);

  // reorder all weights
  SWE->Dx = reorder_2D_fp_arr(SWE->Dx, mapping, SWE->NNodes, SWE->NNbrs_Padded, SWE->NNbrs);
  SWE->Dy = reorder_2D_fp_arr(SWE->Dy, mapping, SWE->NNodes, SWE->NNbrs_Padded, SWE->NNbrs);
  SWE->Dz = reorder_2D_fp_arr(SWE->Dz, mapping, SWE->NNodes, SWE->NNbrs_Padded, SWE->NNbrs);
  SWE->L  = reorder_2D_fp_arr(SWE->L,  mapping, SWE->NNodes, SWE->NNbrs_Padded, SWE->NNbrs);

  // reorder idx
  SWE->idx = reorder_2D_int_arr(SWE->idx, mapping, SWE->NNodes, SWE->NNbrs_Padded, SWE->NNbrs);

  // create inverse mapping: inv_mapping[old_idx] == new_idx
  int* inv_mapping = (int*) malloc(sizeof(int)*SWE->NNodes);
  
  for (int i = 0; i < SWE->NNodes; i++)
    inv_mapping[mapping[i]] = i;

  // account for node id reordering in idx
  for (int i = 0; i < SWE->NNodes; i++){
    for(int j = 0; j < SWE->NNbrs; j++) {
      SWE->idx[i*SWE->NNbrs_Padded+j] = inv_mapping[SWE->idx[i*SWE->NNbrs_Padded+j]];
    }
  }

  // free inverse mapping data
  free(inv_mapping);

}
