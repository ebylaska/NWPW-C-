/*
 * This file generated on line 811 of /home/bylaska/Codes/Scott-DFT/ATLAS/Dev1/..//tune/blas/ger/r1hgen.c
 */
#ifndef ATLAS_DR1KERNELS_H
   #define ATLAS_DR1KERNELS_H

void ATL_dgerk_L0(ATL_CINT, ATL_CINT, const double, const double*, ATL_CINT, const double*, ATL_CINT, double*, ATL_CINT);
void ATL_dgerk_L2(ATL_CINT, ATL_CINT, const double, const double*, ATL_CINT, const double*, ATL_CINT, double*, ATL_CINT);
void ATL_dgerk_L1b(ATL_CINT, ATL_CINT, const double, const double*, ATL_CINT, const double*, ATL_CINT, double*, ATL_CINT);

#define ATL_dgerk_L0_restrict    ATL_dgerk_L0
#define ATL_dgerk_L2_restrict    ATL_dgerk_L2
#define ATL_dgerk_L1             ATL_dgerk_L0
#define ATL_dgerk_L1_restrict    ATL_dgerk_L0
#define ATL_dgerk_L1b_restrict   ATL_dgerk_L1b

#endif /* end guard around atlas_dr1kernels.h */
