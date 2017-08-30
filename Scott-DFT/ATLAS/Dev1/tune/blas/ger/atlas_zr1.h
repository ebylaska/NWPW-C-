#ifndef ATLAS_ZR1_L0_H
#define ATLAS_ZR1_L0_H

#include "atlas_type.h"

#define ATL_r1CacheElts 1536000
#define ATL_r1MU 16
#define ATL_r1NU 2
#define ATL_zgerk_L0 ATL_UGERK
void ATL_zgerk_L0(ATL_CINT, ATL_CINT, const double*, const double*, ATL_CINT, const double*, ATL_CINT, double*, ATL_CINT);
#define ATL_GERK ATL_zgerk_L0
#define ATL_GERKr ATL_zgerk_L0_restrict

#define ATL_GetPartR1(A_, lda_, mb_, nb_) { (mb_) = 307184; (nb_) = ATL_r1NU; }

#endif  /* end protection around header file contents */