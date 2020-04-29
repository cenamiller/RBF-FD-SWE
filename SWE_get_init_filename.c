#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Only need _VERBOSE_ on for debugging
#undef _VERBOSE_

void SWE_get_init_filename(char *Filename, int NNodes, char *testcase){

  char node_str[10];

  if(NNodes >= 1000000000){
    printf("ERROR in SWE_get_init_filename: NNodes to large\n");
    exit(-1);
  }

  strcpy(Filename, "icos");
  snprintf(node_str, 10, "%d", NNodes);
#ifdef _VERBOSE_
  printf("Decimal value of Nnodes = %s\n", node_str);
#endif
  strcat(Filename, node_str);
  strcat(Filename, "_");
  strcat(Filename, testcase);
  strcat(Filename, "_");
  strcat(Filename, "input.bin");

}
