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
#ifndef ATLAS_PTLEVEL2_H
#define ATLAS_PTLEVEL2_H
/*
 * =====================================================================
 * Include files
 * =====================================================================
 */
#include "atlas_enum.h"
/*
 * =====================================================================
 * Prototypes  for single precision real  Level 2  multi-threaded  ATLAS
 * BLAS routines.
 * =====================================================================
 */
void              ATL_sptgbmv
(  const enum ATLAS_TRANS,           const int,       const int,
   const int,       const int,       const float,     const float  *,
   const int,       const float  *,  const int,       const float,
   float  *,        const int );
void              ATL_sptgemv
(  const enum ATLAS_TRANS,           const int,       const int,
   const float,     const float  *,  const int,       const float  *,
   const int,       const float,     float  *,        const int );
void              ATL_sptger
(  const int,       const int,       const float,     const float  *,
   const int,       const float  *,  const int,       float  *,
   const int );
void              ATL_sptsbmv
(  const enum ATLAS_UPLO,            const int,       const int,
   const float,     const float  *,  const int,       const float  *,
   const int,       const float,     float  *,        const int );
void              ATL_sptspmv
(  const enum ATLAS_UPLO,            const int,       const float,
   const float  *,  const float  *,  const int,       const float,
   float  *,        const int );
void              ATL_sptsymv
(  const enum ATLAS_UPLO,            const int,       const float,
   const float  *,  const int,       const float  *,  const int,
   const float,     float  *,        const int );
void              ATL_sptspr
(  const enum ATLAS_UPLO,            const int,       const float,
   const float  *,  const int,       float  * );
void              ATL_sptsyr
(  const enum ATLAS_UPLO,            const int,       const float,
   const float  *,  const int,       float  *,        const int );
void              ATL_sptspr2
(  const enum ATLAS_UPLO,            const int,       const float,
   const float  *,  const int,       const float  *,  const int,
   float  * );
void              ATL_sptsyr2
(  const enum ATLAS_UPLO,            const int,       const float,
   const float  *,  const int,       const float  *,  const int,
   float  *,        const int );
void              ATL_spttbmv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const int,
   const float  *,  const int,       float  *,        const int );
void              ATL_spttpmv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const float  *,
   float  *,        const int );
void              ATL_spttrmv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const float  *,
   const int,       float  *,        const int );
void              ATL_spttbsv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const int,
   const float  *,  const int,       float  *,        const int );
void              ATL_spttpsv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const float  *,
   float  *,        const int );
void              ATL_spttrsv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const float  *,
   const int,       float  *,        const int );
/*
 * =====================================================================
 * Prototypes  for double precision real  Level 2  multi-threaded  ATLAS
 * BLAS routines.
 * =====================================================================
 */
void              ATL_dptgbmv
(  const enum ATLAS_TRANS,           const int,       const int,
   const int,       const int,       const double,    const double *,
   const int,       const double *,  const int,       const double,
   double *,        const int );
void              ATL_dptgemv
(  const enum ATLAS_TRANS,           const int,       const int,
   const double,    const double *,  const int,       const double *,
   const int,       const double,    double *,        const int );
void              ATL_dptger
(  const int,       const int,       const double,    const double *,
   const int,       const double *,  const int,       double *,
   const int );
void              ATL_dptsbmv
(  const enum ATLAS_UPLO,            const int,       const int,
   const double,    const double *,  const int,       const double *,
   const int,       const double,    double *,        const int );
void              ATL_dptspmv
(  const enum ATLAS_UPLO,            const int,       const double,
   const double *,  const double *,  const int,       const double,
   double *,        const int );
void              ATL_dptsymv
(  const enum ATLAS_UPLO,            const int,       const double,
   const double *,  const int,       const double *,  const int,
   const double,    double *,        const int );
void              ATL_dptspr
(  const enum ATLAS_UPLO,            const int,       const double,
   const double *,  const int,       double * );
void              ATL_dptsyr
(  const enum ATLAS_UPLO,            const int,       const double,
   const double *,  const int,       double *,        const int );
void              ATL_dptspr2
(  const enum ATLAS_UPLO,            const int,       const double,
   const double *,  const int,       const double *,  const int,
   double * );
void              ATL_dptsyr2
(  const enum ATLAS_UPLO,            const int,       const double,
   const double *,  const int,       const double *,  const int,
   double *,        const int );
void              ATL_dpttbmv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const int,
   const double *,  const int,       double *,        const int );
void              ATL_dpttpmv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const double *,
   double *,        const int );
void              ATL_dpttrmv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const double *,
   const int,       double *,        const int );
void              ATL_dpttbsv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const int,
   const double *,  const int,       double *,        const int );
void              ATL_dpttpsv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const double *,
   double *,        const int );
void              ATL_dpttrsv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const double *,
   const int,       double *,        const int );
/*
 * =====================================================================
 * Prototypes  for single precision complex Level 2 multi-threaded ATLAS
 * BLAS routines.
 * =====================================================================
 */
void              ATL_cptgbmv
(  const enum ATLAS_TRANS,           const int,       const int,
   const int,       const int,       const float  *,  const float  *,
   const int,       const float  *,  const int,       const float  *,
   float  *,        const int );
void              ATL_cptgemv
(  const enum ATLAS_TRANS,           const int,       const int,
   const float  *,  const float  *,  const int,       const float  *,
   const int,       const float  *,  float  *,        const int );
void              ATL_cptgerc
(  const int,       const int,       const float  *,  const float  *,
   const int,       const float  *,  const int,       float  *,
   const int );
void              ATL_cptgeru
(  const int,       const int,       const float  *,  const float  *,
   const int,       const float  *,  const int,       float  *,
   const int );
void              ATL_cpthbmv
(  const enum ATLAS_UPLO,            const int,       const int,
   const float  *,  const float  *,  const int,       const float  *,
   const int,       const float  *,  float  *,        const int );
void              ATL_cpthpmv
(  const enum ATLAS_UPLO,            const int,       const float  *,
   const float  *,  const float  *,  const int,       const float  *,
   float  *,        const int );
void              ATL_cpthemv
(  const enum ATLAS_UPLO,            const int,       const float  *,
   const float  *,  const int,       const float  *,  const int,
   const float  *,  float  *,        const int );
void              ATL_cpthpr
(  const enum ATLAS_UPLO,            const int,       const float,
   const float  *,  const int,       float  * );
void              ATL_cpther
(  const enum ATLAS_UPLO,            const int,       const float,
   const float  *,  const int,       float  *,        const int );
void              ATL_cpthpr2
(  const enum ATLAS_UPLO,            const int,       const float  *,
   const float  *,  const int,       const float  *,  const int,
   float  * );
void              ATL_cpther2
(  const enum ATLAS_UPLO,            const int,       const float  *,
   const float  *,  const int,       const float  *,  const int,
   float  *,        const int );
void              ATL_cpttbmv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const int,
   const float  *,  const int,       float  *,        const int );
void              ATL_cpttpmv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const float  *,
   float  *,        const int );
void              ATL_cpttrmv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const float  *,
   const int,       float  *,        const int );
void              ATL_cpttbsv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const int,
   const float  *,  const int,       float  *,        const int );
void              ATL_cpttpsv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const float  *,
   float  *,        const int );
void              ATL_cpttrsv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const float  *,
   const int,       float  *,        const int );
/*
 * =====================================================================
 * Prototypes  for double precision complex Level 2 multi-threaded ATLAS
 * BLAS routines.
 * =====================================================================
 */
void              ATL_zptgbmv
(  const enum ATLAS_TRANS,           const int,       const int,
   const int,       const int,       const double *,  const double *,
   const int,       const double *,  const int,       const double *,
   double *,        const int );
void              ATL_zptgemv
(  const enum ATLAS_TRANS,           const int,       const int,
   const double *,  const double *,  const int,       const double *,
   const int,       const double *,  double *,        const int );
void              ATL_zptgerc
(  const int,       const int,       const double *,  const double *,
   const int,       const double *,  const int,       double *,
   const int );
void              ATL_zptgeru
(  const int,       const int,       const double *,  const double *,
   const int,       const double *,  const int,       double *,
   const int );
void              ATL_zpthbmv
(  const enum ATLAS_UPLO,            const int,       const int,
   const double *,  const double *,  const int,       const double *,
   const int,       const double *,  double *,        const int );
void              ATL_zpthpmv
(  const enum ATLAS_UPLO,            const int,       const double *,
   const double *,  const double *,  const int,       const double *,
   double *,        const int );
void              ATL_zpthemv
(  const enum ATLAS_UPLO,            const int,       const double *,
   const double *,  const int,       const double *,  const int,
   const double *,  double *,        const int );
void              ATL_zpthpr
(  const enum ATLAS_UPLO,            const int,       const double,
   const double *,  const int,       double * );
void              ATL_zpther
(  const enum ATLAS_UPLO,            const int,       const double,
   const double *,  const int,       double *,        const int );
void              ATL_zpthpr2
(  const enum ATLAS_UPLO,            const int,       const double *,
   const double *,  const int,       const double *,  const int,
   double * );
void              ATL_zpther2
(  const enum ATLAS_UPLO,            const int,       const double *,
   const double *,  const int,       const double *,  const int,
   double *,        const int );
void              ATL_zpttbmv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const int,
   const double *,  const int,       double *,        const int );
void              ATL_zpttpmv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const double *,
   double *,        const int );
void              ATL_zpttrmv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const double *,
   const int,       double *,        const int );
void              ATL_zpttbsv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const int,
   const double *,  const int,       double *,        const int );
void              ATL_zpttpsv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const double *,
   double *,        const int );
void              ATL_zpttrsv
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const double *,
   const int,       double *,        const int );

#endif
/*
 * End of atlas_ptlevel2.h
 */
