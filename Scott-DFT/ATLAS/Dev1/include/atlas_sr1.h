#ifndef ATLAS_SR1_L0_H
#define ATLAS_SR1_L0_H

#include "atlas_type.h"

#define ATL_r1CacheElts 0
#define ATL_r1MU 16
#define ATL_r1NU 1
#define ATL_r1NOBLOCK
#include "atlas_sr1kernels.h"
#define ATL_GERK ATL_sgerk_L0
#define ATL_GERKr ATL_sgerk_L0_restrict

#define ATL_GetPartR1(A_, lda_, mb_, nb_) { (mb_) = 0; (nb_) = ATL_r1NU; }

#endif  /* end protection around header file contents */
