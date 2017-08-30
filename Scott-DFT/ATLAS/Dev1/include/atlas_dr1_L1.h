#ifndef ATLAS_DR1_L1_H
#define ATLAS_DR1_L1_H

#include "atlas_type.h"

#define ATL_r1CacheElts 10240
#define ATL_r1MU 16
#define ATL_r1NU 1
#include "atlas_dr1kernels.h"
#define ATL_GERK ATL_dgerk_L1
#define ATL_GERKr ATL_dgerk_L1_restrict

#define ATL_GetPartR1(A_, lda_, mb_, nb_) { (mb_) = 3408; (nb_) = ATL_r1NU; }

#endif  /* end protection around header file contents */
