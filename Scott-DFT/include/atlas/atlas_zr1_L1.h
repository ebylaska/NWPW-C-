#ifndef ATLAS_ZR1_L1_H
#define ATLAS_ZR1_L1_H

#include "atlas_type.h"

#define ATL_r1CacheElts 3072
#define ATL_r1MU 2
#define ATL_r1NU 1
#include "atlas_zr1kernels.h"
#define ATL_GERK ATL_zgerk_L1
#define ATL_GERKr ATL_zgerk_L1_restrict

#define ATL_GetPartR1(A_, lda_, mb_, nb_) { (mb_) = 1022; (nb_) = ATL_r1NU; }

#endif  /* end protection around header file contents */
