/*
 * This file generated on line 811 of /home/bylaska/Codes/Scott-DFT/ATLAS/Dev1/..//tune/blas/ger/r1hgen.c
 */
#ifndef ATLAS_ZR1KERNELS_H
   #define ATLAS_ZR1KERNELS_H

void ATL_zgerk_L0(ATL_CINT, ATL_CINT, const double*, const double*, ATL_CINT, const double*, ATL_CINT, double*, ATL_CINT);
void ATL_zgerk_L2(ATL_CINT, ATL_CINT, const double*, const double*, ATL_CINT, const double*, ATL_CINT, double*, ATL_CINT);
void ATL_zgerk_L1(ATL_CINT, ATL_CINT, const double*, const double*, ATL_CINT, const double*, ATL_CINT, double*, ATL_CINT);

#define ATL_zgerk_L0_restrict    ATL_zgerk_L0
#define ATL_zgerk_L2_restrict    ATL_zgerk_L2
#define ATL_zgerk_L1_restrict    ATL_zgerk_L1
#define ATL_zgerk_L1b            ATL_zgerk_L0
#define ATL_zgerk_L1b_restrict   ATL_zgerk_L0

#endif /* end guard around atlas_zr1kernels.h */
