#ifndef ATLAS_CR1_L2_H
#define ATLAS_CR1_L2_H

#include "atlas_type.h"

#define ATL_r1CacheElts 2129
#define ATL_r1MU 16
#define ATL_r1NU 1
#include "atlas_cr1kernels.h"
#define ATL_GERK ATL_cgerk_L2
#define ATL_GERKr ATL_cgerk_L2_restrict

#define ATL_GetPartR1(A_, lda_, mb_, nb_) { (mb_) = 704; (nb_) = ATL_r1NU; }

#endif  /* end protection around header file contents */
