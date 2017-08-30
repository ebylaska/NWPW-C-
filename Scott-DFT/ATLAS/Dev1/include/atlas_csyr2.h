/*
 * This file generated on line 603 of /home/bylaska/Codes/Scott-DFT/ATLAS/Dev1/..//tune/blas/ger/r1hgen.c
 */
#ifndef ATLAS_CSYR2_H
   #define ATLAS_CSYR2_H

#include "atlas_cr1kernels.h"
#define ATL_s2CacheElts 32768
#define ATL_s2MU 16
#define ATL_s2NU 1
#define ATL_R1OC ATL_cgerk_L0
#define ATL_R1OCr ATL_cgerk_L0
#define ATL_R1IC ATL_cgerk_L2
#define ATL_R1ICr ATL_cgerk_L2

#define ATL_GetPartS2(A_, lda_, mb_, nb_) { (mb_) = 8176; (nb_) = ATL_s2NU; }

#define ATL_HER2U_nu(A_, lda_, x_, y_, xt_, yt_) \
{ \
   TYPE *aa=(A_); \
   ATL_CINT lda0_ = 0; \
   const TYPE x0r=*(x_), x0i=(x_)[1]; \
   const TYPE xt0r=*(xt_), xt0i=(xt_)[1]; \
   const TYPE y0r=*(y_), y0i=(y_)[1]; \
   const TYPE yt0r=*(yt_), yt0i=(yt_)[1]; \
   aa[lda0_+0] += x0r*yt0r-x0i*yt0i + y0r*xt0r-y0i*xt0i; \
   aa[lda0_+1] = 0.0; \
}
#define ATL_HER2L_nu(A_, lda_, x_, y_, xt_, yt_) \
{ \
   TYPE *aa=(A_); \
   ATL_CINT lda0_ = 0; \
   const TYPE x0r=*(x_), x0i=(x_)[1]; \
   const TYPE xt0r=*(xt_), xt0i=(xt_)[1]; \
   const TYPE y0r=*(y_), y0i=(y_)[1]; \
   const TYPE yt0r=*(yt_), yt0i=(yt_)[1]; \
   aa[lda0_+0] += x0r*yt0r-x0i*yt0i + y0r*xt0r-y0i*xt0i; \
   aa[lda0_+1] = 0.0; \
}

#endif
