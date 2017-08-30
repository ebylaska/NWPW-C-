#ifndef ATLAS_DR1_L2_H
#define ATLAS_DR1_L2_H

#include "atlas_type.h"

#define ATL_r1CacheElts 2048
#define ATL_r1MU 4
#define ATL_r1NU 4
#include "atlas_dr1kernels.h"
#define ATL_GERK ATL_dgerk_L2
#define ATL_GERKr ATL_dgerk_L2_restrict

#define ATL_GetPartR1(A_, lda_, mb_, nb_) { (mb_) = 224; (nb_) = ATL_r1NU; }

#endif  /* end protection around header file contents */
