#ifndef ATLAS_SR1_L2_H
#define ATLAS_SR1_L2_H

#include "atlas_type.h"

#define ATL_r1CacheElts 4915
#define ATL_r1MU 4
#define ATL_r1NU 4
#include "atlas_sr1kernels.h"
#define ATL_GERK ATL_sgerk_L2
#define ATL_GERKr ATL_sgerk_L2_restrict

#define ATL_GetPartR1(A_, lda_, mb_, nb_) { (mb_) = 544; (nb_) = ATL_r1NU; }

#endif  /* end protection around header file contents */
