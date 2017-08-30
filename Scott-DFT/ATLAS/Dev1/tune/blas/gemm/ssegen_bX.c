#define ATL_INT int
#include <stdio.h>
#include <stdint.h>
#include <pmmintrin.h>

#define I_UNROLL 4
#define J_UNROLL 3
/* Is prefetched data written or just read? */
#define PF_READONLY 0
#define PF_READWRITE 1
#define PF_NO_REUSE 0

/* Default temporality of cache prefetch (1-3) */
#define PF_DEF 1
#define CACHE_LINE_SIZE 64
#define MMCAST( a ) (a)
#define MMCASTintrin( a ) (a)
#define TYPE float
void ATL_USERMM( const ATL_INT M, const ATL_INT N, const ATL_INT K,
                 const TYPE alpha, const TYPE *A, const ATL_INT lda,
                 const TYPE *B, const ATL_INT ldb,
                 const TYPE beta, TYPE *C, const ATL_INT ldc )
{
   register ATL_INT i, j, k;
   
   /* Vector registers to hold the elements of C */
   __m128 c0_0, c0_1, c0_2, c0_3;
   __m128 c1_0, c1_1, c1_2, c1_3;
   __m128 c2_0, c2_1, c2_2, c2_3;
   /* Vector register to hold C*beta */
   __m128 bc0_0, bc0_1, bc0_2, bc0_3;
   __m128 bc1_0, bc1_1, bc1_2, bc1_3;
   __m128 bc2_0, bc2_1, bc2_2, bc2_3;
   /* Temporary vector registers for use in inner loop */
   __m128 temp; 
   __m128 temp0;  
   __m128 temp1;  
   __m128 temp2;  
   __m128 temp3;  
   const __m128 betaV = _mm_set1_ps( beta ); 
   /* Pointer adjustments */  
   register const ATL_INT ldc_bytes = 2*ldc;
   
   register TYPE const *B0_off = B;
      
   
   /* Unroll A */
   __m128 A0, a0, A1, a1, A2, a2, A3, a3;
   /* Unroll B */
   __m128 B0, B1, B2;
   

   register const ATL_INT unroll_a = I_UNROLL*KB;
   register TYPE* cPtr = C;
   

   /* =======================================
    * Begin generated inner loops for case Non aligned
    * ======================================= */
   for( j=-72; j!=0; j+=J_UNROLL) 
   {
      register TYPE const *A0_off = A; 
      
      register TYPE *cPtrI0 = cPtr;
      register TYPE *cPtrI1 = cPtrI0 + ldc_bytes;
      register TYPE *cPtrI2 = cPtrI1 + ldc_bytes;
      

      for( i=-72; i != 0; i+= I_UNROLL )
      {
         /* K_Unrolling0 */
         A0 = _mm_load_ps( MMCAST(A0_off) );
         A1 = _mm_load_ps( MMCAST(A0_off + 72) );
         A2 = _mm_load_ps( MMCAST(A0_off + 144) );
         A3 = _mm_load_ps( MMCAST(A0_off + 216) );
         B0 = _mm_load_ps( MMCAST(B0_off) );
         c0_0 = B0;
         c0_0 = _mm_mul_ps( A0, c0_0 );
         c0_1 = B0;
         c0_1 = _mm_mul_ps( A1, c0_1 );
         c0_2 = B0;
         c0_2 = _mm_mul_ps( A2, c0_2 );
         c0_3 = B0;
         c0_3 = _mm_mul_ps( A3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 72) );
         c1_0 = B1;
         c1_0 = _mm_mul_ps( A0, c1_0 );
         c1_1 = B1;
         c1_1 = _mm_mul_ps( A1, c1_1 );
         c1_2 = B1;
         c1_2 = _mm_mul_ps( A2, c1_2 );
         c1_3 = B1;
         c1_3 = _mm_mul_ps( A3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 144) );
         c2_0 = B2;
         c2_0 = _mm_mul_ps( A0, c2_0 );
         c2_1 = B2;
         c2_1 = _mm_mul_ps( A1, c2_1 );
         c2_2 = B2;
         c2_2 = _mm_mul_ps( A2, c2_2 );
         c2_3 = B2;
         c2_3 = _mm_mul_ps( A3, c2_3 );
         
         /* K_Unrolling: 4 */
         A0 = _mm_load_ps( MMCAST(A0_off + 4) );
         A1 = _mm_load_ps( MMCAST(A0_off + 76) );
         A2 = _mm_load_ps( MMCAST(A0_off + 148) );
         A3 = _mm_load_ps( MMCAST(A0_off + 220) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 4) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 76) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 148) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* K_Unrolling: 8 */
         A0 = _mm_load_ps( MMCAST(A0_off + 8) );
         A1 = _mm_load_ps( MMCAST(A0_off + 80) );
         A2 = _mm_load_ps( MMCAST(A0_off + 152) );
         A3 = _mm_load_ps( MMCAST(A0_off + 224) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 8) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 80) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 152) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* K_Unrolling: 12 */
         A0 = _mm_load_ps( MMCAST(A0_off + 12) );
         A1 = _mm_load_ps( MMCAST(A0_off + 84) );
         A2 = _mm_load_ps( MMCAST(A0_off + 156) );
         A3 = _mm_load_ps( MMCAST(A0_off + 228) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 12) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 84) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 156) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* K_Unrolling: 16 */
         A0 = _mm_load_ps( MMCAST(A0_off + 16) );
         A1 = _mm_load_ps( MMCAST(A0_off + 88) );
         A2 = _mm_load_ps( MMCAST(A0_off + 160) );
         A3 = _mm_load_ps( MMCAST(A0_off + 232) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 16) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 88) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 160) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* K_Unrolling: 20 */
         A0 = _mm_load_ps( MMCAST(A0_off + 20) );
         A1 = _mm_load_ps( MMCAST(A0_off + 92) );
         A2 = _mm_load_ps( MMCAST(A0_off + 164) );
         A3 = _mm_load_ps( MMCAST(A0_off + 236) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 20) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 92) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 164) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* K_Unrolling: 24 */
         A0 = _mm_load_ps( MMCAST(A0_off + 24) );
         A1 = _mm_load_ps( MMCAST(A0_off + 96) );
         A2 = _mm_load_ps( MMCAST(A0_off + 168) );
         A3 = _mm_load_ps( MMCAST(A0_off + 240) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 24) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 96) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 168) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* K_Unrolling: 28 */
         A0 = _mm_load_ps( MMCAST(A0_off + 28) );
         A1 = _mm_load_ps( MMCAST(A0_off + 100) );
         A2 = _mm_load_ps( MMCAST(A0_off + 172) );
         A3 = _mm_load_ps( MMCAST(A0_off + 244) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 28) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 100) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 172) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* K_Unrolling: 32 */
         A0 = _mm_load_ps( MMCAST(A0_off + 32) );
         A1 = _mm_load_ps( MMCAST(A0_off + 104) );
         A2 = _mm_load_ps( MMCAST(A0_off + 176) );
         A3 = _mm_load_ps( MMCAST(A0_off + 248) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 32) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 104) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 176) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* K_Unrolling: 36 */
         A0 = _mm_load_ps( MMCAST(A0_off + 36) );
         A1 = _mm_load_ps( MMCAST(A0_off + 108) );
         A2 = _mm_load_ps( MMCAST(A0_off + 180) );
         A3 = _mm_load_ps( MMCAST(A0_off + 252) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 36) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 108) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 180) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* K_Unrolling: 40 */
         A0 = _mm_load_ps( MMCAST(A0_off + 40) );
         A1 = _mm_load_ps( MMCAST(A0_off + 112) );
         A2 = _mm_load_ps( MMCAST(A0_off + 184) );
         A3 = _mm_load_ps( MMCAST(A0_off + 256) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 40) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 112) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 184) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* K_Unrolling: 44 */
         A0 = _mm_load_ps( MMCAST(A0_off + 44) );
         A1 = _mm_load_ps( MMCAST(A0_off + 116) );
         A2 = _mm_load_ps( MMCAST(A0_off + 188) );
         A3 = _mm_load_ps( MMCAST(A0_off + 260) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 44) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 116) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 188) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* K_Unrolling: 48 */
         A0 = _mm_load_ps( MMCAST(A0_off + 48) );
         A1 = _mm_load_ps( MMCAST(A0_off + 120) );
         A2 = _mm_load_ps( MMCAST(A0_off + 192) );
         A3 = _mm_load_ps( MMCAST(A0_off + 264) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 48) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 120) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 192) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* K_Unrolling: 52 */
         A0 = _mm_load_ps( MMCAST(A0_off + 52) );
         A1 = _mm_load_ps( MMCAST(A0_off + 124) );
         A2 = _mm_load_ps( MMCAST(A0_off + 196) );
         A3 = _mm_load_ps( MMCAST(A0_off + 268) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 52) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 124) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 196) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* K_Unrolling: 56 */
         A0 = _mm_load_ps( MMCAST(A0_off + 56) );
         A1 = _mm_load_ps( MMCAST(A0_off + 128) );
         A2 = _mm_load_ps( MMCAST(A0_off + 200) );
         A3 = _mm_load_ps( MMCAST(A0_off + 272) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 56) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 128) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 200) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* K_Unrolling: 60 */
         A0 = _mm_load_ps( MMCAST(A0_off + 60) );
         A1 = _mm_load_ps( MMCAST(A0_off + 132) );
         A2 = _mm_load_ps( MMCAST(A0_off + 204) );
         A3 = _mm_load_ps( MMCAST(A0_off + 276) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 60) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 132) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 204) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* K_Unrolling: 64 */
         A0 = _mm_load_ps( MMCAST(A0_off + 64) );
         A1 = _mm_load_ps( MMCAST(A0_off + 136) );
         A2 = _mm_load_ps( MMCAST(A0_off + 208) );
         A3 = _mm_load_ps( MMCAST(A0_off + 280) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 64) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 136) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 208) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* K_Unrolling: 68 */
         A0 = _mm_load_ps( MMCAST(A0_off + 68) );
         A1 = _mm_load_ps( MMCAST(A0_off + 140) );
         A2 = _mm_load_ps( MMCAST(A0_off + 212) );
         A3 = _mm_load_ps( MMCAST(A0_off + 284) );
         
         B0 = _mm_load_ps( MMCAST(B0_off + 68) );
         a0 = A0;
         a0 = _mm_mul_ps( B0, a0 );
         c0_0 = _mm_add_ps( a0, c0_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B0, a1 );
         c0_1 = _mm_add_ps( a1, c0_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B0, a2 );
         c0_2 = _mm_add_ps( a2, c0_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B0, a3 );
         c0_3 = _mm_add_ps( a3, c0_3 );
         
         B1 = _mm_load_ps( MMCAST(B0_off + 140) );
         a0 = A0;
         a0 = _mm_mul_ps( B1, a0 );
         c1_0 = _mm_add_ps( a0, c1_0 );
         a1 = A1;
         a1 = _mm_mul_ps( B1, a1 );
         c1_1 = _mm_add_ps( a1, c1_1 );
         a2 = A2;
         a2 = _mm_mul_ps( B1, a2 );
         c1_2 = _mm_add_ps( a2, c1_2 );
         a3 = A3;
         a3 = _mm_mul_ps( B1, a3 );
         c1_3 = _mm_add_ps( a3, c1_3 );
         
         B2 = _mm_load_ps( MMCAST(B0_off + 212) );
         A0 = _mm_mul_ps( B2, A0 );
         c2_0 = _mm_add_ps( A0, c2_0 );
         A1 = _mm_mul_ps( B2, A1 );
         c2_1 = _mm_add_ps( A1, c2_1 );
         A2 = _mm_mul_ps( B2, A2 );
         c2_2 = _mm_add_ps( A2, c2_2 );
         A3 = _mm_mul_ps( B2, A3 );
         c2_3 = _mm_add_ps( A3, c2_3 );
         /* Single Scalar Compression */
         c0_0 = _mm_hadd_ps( c0_0, c0_1 );
         c0_2 = _mm_hadd_ps( c0_2, c0_3 );
         c0_0 = _mm_hadd_ps( c0_0, c0_2 );
         
         c1_0 = _mm_hadd_ps( c1_0, c1_1 );
         c1_2 = _mm_hadd_ps( c1_2, c1_3 );
         c1_0 = _mm_hadd_ps( c1_0, c1_2 );
         
         c2_0 = _mm_hadd_ps( c2_0, c2_1 );
         c2_2 = _mm_hadd_ps( c2_2, c2_3 );
         c2_0 = _mm_hadd_ps( c2_0, c2_2 );
         
         /* Applying Beta */
            /* Apply Beta Factor */
            /* Load C from memory */
            temp0 = _mm_load_ss( cPtrI0 + 0 );
            temp1 = _mm_load_ss( cPtrI0 + 2 );
            temp2 = _mm_load_ss( cPtrI0 + 4 );
            temp3 = _mm_load_ss( cPtrI0 + 6 );
            temp0 = _mm_unpacklo_ps( temp0, temp1 );
            temp2 = _mm_unpacklo_ps( temp2, temp3 );
            bc0_0 = _mm_movelh_ps( temp0, temp2 );
            bc0_0 = _mm_mul_ps( betaV, bc0_0 );
            /* Load C from memory */
            temp0 = _mm_load_ss( cPtrI1 + 0 );
            temp1 = _mm_load_ss( cPtrI1 + 2 );
            temp2 = _mm_load_ss( cPtrI1 + 4 );
            temp3 = _mm_load_ss( cPtrI1 + 6 );
            temp0 = _mm_unpacklo_ps( temp0, temp1 );
            temp2 = _mm_unpacklo_ps( temp2, temp3 );
            bc1_0 = _mm_movelh_ps( temp0, temp2 );
            bc1_0 = _mm_mul_ps( betaV, bc1_0 );
            /* Load C from memory */
            temp0 = _mm_load_ss( cPtrI2 + 0 );
            temp1 = _mm_load_ss( cPtrI2 + 2 );
            temp2 = _mm_load_ss( cPtrI2 + 4 );
            temp3 = _mm_load_ss( cPtrI2 + 6 );
            temp0 = _mm_unpacklo_ps( temp0, temp1 );
            temp2 = _mm_unpacklo_ps( temp2, temp3 );
            bc2_0 = _mm_movelh_ps( temp0, temp2 );
            bc2_0 = _mm_mul_ps( betaV, bc2_0 );
            /* C = (beta*C) + (matrix multiply) */
            c0_0 = _mm_add_ps( bc0_0, c0_0 );
            c1_0 = _mm_add_ps( bc1_0, c1_0 );
            c2_0 = _mm_add_ps( bc2_0, c2_0 );
         /* Move pointers to next iteration */  
         A0_off += unroll_a;
         
         /* Store results back to memory  */
         temp = c0_0;
         _mm_store_ss( cPtrI0+0, temp );
         temp = _mm_shuffle_ps( c0_0, c0_0,_MM_SHUFFLE(1,1,1,1));
         _mm_store_ss( cPtrI0+2, temp );
         temp = _mm_shuffle_ps( c0_0, c0_0,_MM_SHUFFLE(2,2,2,2));
         _mm_store_ss( cPtrI0+4, temp );
         temp = _mm_shuffle_ps( c0_0, c0_0,_MM_SHUFFLE(3,3,3,3));
         _mm_store_ss( cPtrI0+6, temp );
         temp = c1_0;
         _mm_store_ss( cPtrI1+0, temp );
         temp = _mm_shuffle_ps( c1_0, c1_0,_MM_SHUFFLE(1,1,1,1));
         _mm_store_ss( cPtrI1+2, temp );
         temp = _mm_shuffle_ps( c1_0, c1_0,_MM_SHUFFLE(2,2,2,2));
         _mm_store_ss( cPtrI1+4, temp );
         temp = _mm_shuffle_ps( c1_0, c1_0,_MM_SHUFFLE(3,3,3,3));
         _mm_store_ss( cPtrI1+6, temp );
         temp = c2_0;
         _mm_store_ss( cPtrI2+0, temp );
         temp = _mm_shuffle_ps( c2_0, c2_0,_MM_SHUFFLE(1,1,1,1));
         _mm_store_ss( cPtrI2+2, temp );
         temp = _mm_shuffle_ps( c2_0, c2_0,_MM_SHUFFLE(2,2,2,2));
         _mm_store_ss( cPtrI2+4, temp );
         temp = _mm_shuffle_ps( c2_0, c2_0,_MM_SHUFFLE(3,3,3,3));
         _mm_store_ss( cPtrI2+6, temp );
         cPtrI0 += 2*I_UNROLL;
         cPtrI1 += 2*I_UNROLL;
         cPtrI2 += 2*I_UNROLL;
         

      } /* End i/MB loop */

      B0_off += J_UNROLL*KB;
      cPtr += J_UNROLL*ldc_bytes;
   } /* End j/NB loop */
   /* End of generated inner loops */
}
