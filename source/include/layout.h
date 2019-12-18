#ifndef SWE_LAYOUT_H
#define SWE_LAYOUT_H

#include <swe_config.h>

GSMD_struct convert2_SFDL(GSMD_struct GSMD);

// pads all GSMD data to padded_Nnodes which aligns with PADDED_NUM
// fills new nodes with the same data as the last node -> simply recalculates these values and doesn't cause any parallelization/distribution issues
GSMD_struct pad_GSMD_data(GSMD_struct GSMD);

#endif
