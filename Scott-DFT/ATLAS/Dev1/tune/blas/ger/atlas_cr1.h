#ifndef ATLAS_CR1_L0_H
#define ATLAS_CR1_L0_H

#include "atlas_type.h"

#define ATL_r1CacheElts 3072000
#define ATL_r1MU 16
#define ATL_r1NU 2
#define ATL_cgerk_L0 ATL_UGERK
void ATL_cgerk_L0(ATL_CINT, ATL_CINT, const float*, const float*, ATL_CINT, const float*, ATL_CINT, float*, ATL_CINT);
#define ATL_GERK ATL_cgerk_L0
#define ATL_GERKr ATL_cgerk_L0_restrict

#define ATL_GetPartR1(A_, lda_, mb_, nb_) { (mb_) = 614384; (nb_) = ATL_r1NU; }

#endif  /* end protection around header file contents */
