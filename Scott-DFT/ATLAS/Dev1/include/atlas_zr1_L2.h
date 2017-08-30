#ifndef ATLAS_ZR1_L2_H
#define ATLAS_ZR1_L2_H

#include "atlas_type.h"

#define ATL_r1CacheElts 1024
#define ATL_r1MU 1
#define ATL_r1NU 4
#include "atlas_zr1kernels.h"
#define ATL_GERK ATL_zgerk_L2
#define ATL_GERKr ATL_zgerk_L2_restrict

#define ATL_GetPartR1(A_, lda_, mb_, nb_) { (mb_) = 112; (nb_) = ATL_r1NU; }

#endif  /* end protection around header file contents */
