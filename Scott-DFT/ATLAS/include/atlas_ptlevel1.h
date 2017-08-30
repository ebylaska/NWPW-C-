/* ---------------------------------------------------------------------
 *
 * -- Automatically Tuned Linear Algebra Software (ATLAS)
 *    (C) Copyright 2000 All Rights Reserved
 *
 * -- ATLAS routine -- Version 3.2 -- December 25, 2000
 *
 * Author         : Antoine P. Petitet
 * Originally developed at the University of Tennessee,
 * Innovative Computing Laboratory, Knoxville TN, 37996-1301, USA.
 *
 * ---------------------------------------------------------------------
 *
 * -- Copyright notice and Licensing terms:
 *
 *  Redistribution  and  use in  source and binary forms, with or without
 *  modification, are  permitted provided  that the following  conditions
 *  are met:
 *
 * 1. Redistributions  of  source  code  must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce  the above copyright
 *    notice,  this list of conditions, and the  following disclaimer in
 *    the documentation and/or other materials provided with the distri-
 *    bution.
 * 3. The name of the University,  the ATLAS group,  or the names of its
 *    contributors  may not be used to endorse or promote products deri-
 *    ved from this software without specific written permission.
 *
 * -- Disclaimer:
 *
 * THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,  INDIRECT, INCIDENTAL, SPE-
 * CIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO,  PROCUREMENT  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEO-
 * RY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  (IN-
 * CLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ---------------------------------------------------------------------
 */
#ifndef ATL_PTLEVEL1_H
#define ATL_PTLEVEL1_H
/*
 * =====================================================================
 * Prototypes  for single precision real  Level 1  multi-threaded  ATLAS
 * BLAS routines.
 * =====================================================================
 */
void              ATL_sptrotg
(  float  *,        float  *,        float  *,        float  * );
void              ATL_sptfrotmg
(  float  *,        float  *,        float  *,        const float,
   float  * );
float             ATL_sptnrm2
(  const int,       const float  *,  const int );
float             ATL_sptasum
(  const int,       const float  *,  const int );
int               ATL_isptamax
(  const int,       const float  *,  const int );
void              ATL_sptscal
(  const int,       const float,     float  *,        const int );
void              ATL_sptaxpy
(  const int,       const float,     const float  *,  const int,
   float  *,        const int );
void              ATL_sptcopy
(  const int,       const float  *,  const int,       float  *,
   const int );
void              ATL_sptswap
(  const int,       float  *,        const int,       float  *,
   const int );
void              ATL_sptrot
(  const int,       float  *,        const int,       float  *,
   const int,       const float,     const float  );
void              ATL_sptrotm
(  const int,       float  *,        const int,       float  *,
   const int,       const float  * );
float             ATL_sptdot
(  const int,       const float  *,  const int,       const float  *,
   const int );
float             ATL_sdsptdot
(  const int,       const float,     const float  *,  const int,
   const float  *,  const int );
/*
 * =====================================================================
 * Prototypes  for double precision real  Level 1  multi-threaded  ATLAS
 * BLAS routines.
 * =====================================================================
 */
void              ATL_dptrotg
(  double *,        double *,        double *,        double * );
void              ATL_dptfrotmg
(  double *,        double *,        double *,        const double,
   double * );
double            ATL_dptnrm2
(  const int,       const double *,  const int );
double            ATL_dptasum
(  const int,       const double *,  const int );
int               ATL_idptamax
(  const int,       const double *,  const int );
void              ATL_dptscal
(  const int,       const double,    double *,        const int );
void              ATL_dptaxpy
(  const int,       const double,    const double *,  const int,
   double *,        const int );
void              ATL_dptcopy
(  const int,       const double *,  const int,       double *,
   const int );
void              ATL_dptswap
(  const int,       double *,        const int,       double *,
   const int );
void              ATL_dptrot
(  const int,       double *,        const int,       double *,
   const int,       const double,    const double );
void              ATL_dptrotm
(  const int,       double *,        const int,       double *,
   const int,       const double * );
double            ATL_dptdot
(  const int,       const double *,  const int,       const double *,
   const int );
double            ATL_dsptdot
(  const int,       const float  *,  const int,       const float  *,
   const int );
/*
 * =====================================================================
 * Prototypes  for single precision complex Level 1 multi-threaded ATLAS
 * BLAS routines.
 * =====================================================================
 */
void              ATL_cptrotg
(  float  *,        const float  *,  float  *,        float  * );
float             ATL_scptnrm2
(  const int,       const float  *,  const int );
float             ATL_scptasum
(  const int,       const float  *,  const int );
int               ATL_icptamax
(  const int,       const float  *,  const int );
void              ATL_cptscal
(  const int,       const float  *,  float  *,        const int );
void              ATL_csptscal
(  const int,       const float,     float  *,        const int );
void              ATL_cptaxpy
(  const int,       const float  *,  const float  *,  const int,
   float  *,        const int );
void              ATL_cptcopy
(  const int,       const float  *,  const int,       float  *,
   const int );
void              ATL_cptswap
(  const int,       float  *,        const int,       float  *,
   const int );
void              ATL_csptrot
(  const int,       float  *,        const int,       float  *,
   const int,       const float,     const float  );
void              ATL_cptdotu_sub
(  const int,       const float  *,  const int,       const float  *,
   const int,       float * );
void              ATL_cptdotc_sub
(  const int,       const float  *,  const int,       const float  *,
   const int,       float * );
/*
 * =====================================================================
 * Prototypes  for single precision complex Level 1 multi-threaded ATLAS
 * BLAS routines.
 * =====================================================================
 */
void              ATL_zptfrotg
(  double *,        const double *,  double *,        double * );
double            ATL_dzptnrm2
(  const int,       const double *,  const int );
double            ATL_dzptasum
(  const int,       const double *,  const int );
int               ATL_izptamax
(  const int,       const double *,  const int );
void              ATL_zptscal
(  const int,       const double *,  double *,        const int );
void              ATL_zdptscal
(  const int,       const double,    double *,        const int );
void              ATL_zptaxpy
(  const int,       const double *,  const double *,  const int,
   double *,        const int );
void              ATL_zptcopy
(  const int,       const double *,  const int,       double *,
   const int );
void              ATL_zptswap
(  const int,       double *,        const int,       double *,
   const int );
void              ATL_zdptrot
(  const int,       double *,        const int,       double *,
   const int,       const double,    const double );
void              ATL_zptdotu_sub
(  const int,       const double *,  const int,       const double *,
   const int,       double * );
void              ATL_zptdotc_sub
(  const int,       const double *,  const int,       const double *,
   const int,       double * );

#endif
/*
 * End of atlas_ptlevel1.h
 */
