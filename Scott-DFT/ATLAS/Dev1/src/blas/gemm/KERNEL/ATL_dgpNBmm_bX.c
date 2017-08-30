#ifdef ATL_UCLEANN
#define ATL_dJIK48x0x48TN48x48x0_a1_bX ATL_dgpNBmm_bX
#else
#define ATL_dJIK48x0x48TN48x48x0_a1_bX ATL_dpNBmm_bX
#endif

#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
static void ATL_dJIK48x0x48TN4x1x48_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=48, NB=0, KB=48, 
 * lda=48, ldb=48, ldc=0, mu=4, nu=1, ku=48, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   #define Nb N
   const double *stM = A + 2304;
   const double *stN = B + (((Nb) << 5)+((Nb) << 4));
   #define incAk 48
   const int incAm = 144, incAn = -2304;
   #define incBk 48
   const int incBm = -48, incBn = 48;
   #define incCm 4
   const int incCn = (ldc) - 48;
   double *pC0=C;
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0, rA1, rA2, rA3;
   register double rB0;
   register double rC0_0, rC1_0, rC2_0, rC3_0;
   do /* N-loop */
   {
      do /* M-loop */
      {
         rA0 = beta;
         rC0_0 = *pC0;
         rC0_0 *= rA0;
         rC1_0 = pC0[1];
         rC1_0 *= rA0;
         rC2_0 = pC0[2];
         rC2_0 *= rA0;
         rC3_0 = pC0[3];
         rC3_0 *= rA0;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[48];
         rA2 = pA0[96];
         rA3 = pA0[144];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[1];
         rB0 = pB0[1];
         rA1 = pA0[49];
         rA2 = pA0[97];
         rA3 = pA0[145];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[2];
         rB0 = pB0[2];
         rA1 = pA0[50];
         rA2 = pA0[98];
         rA3 = pA0[146];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[3];
         rB0 = pB0[3];
         rA1 = pA0[51];
         rA2 = pA0[99];
         rA3 = pA0[147];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[4];
         rB0 = pB0[4];
         rA1 = pA0[52];
         rA2 = pA0[100];
         rA3 = pA0[148];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[5];
         rB0 = pB0[5];
         rA1 = pA0[53];
         rA2 = pA0[101];
         rA3 = pA0[149];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[6];
         rB0 = pB0[6];
         rA1 = pA0[54];
         rA2 = pA0[102];
         rA3 = pA0[150];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[7];
         rB0 = pB0[7];
         rA1 = pA0[55];
         rA2 = pA0[103];
         rA3 = pA0[151];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[8];
         rB0 = pB0[8];
         rA1 = pA0[56];
         rA2 = pA0[104];
         rA3 = pA0[152];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[9];
         rB0 = pB0[9];
         rA1 = pA0[57];
         rA2 = pA0[105];
         rA3 = pA0[153];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[10];
         rB0 = pB0[10];
         rA1 = pA0[58];
         rA2 = pA0[106];
         rA3 = pA0[154];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[11];
         rB0 = pB0[11];
         rA1 = pA0[59];
         rA2 = pA0[107];
         rA3 = pA0[155];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[12];
         rB0 = pB0[12];
         rA1 = pA0[60];
         rA2 = pA0[108];
         rA3 = pA0[156];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[13];
         rB0 = pB0[13];
         rA1 = pA0[61];
         rA2 = pA0[109];
         rA3 = pA0[157];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[14];
         rB0 = pB0[14];
         rA1 = pA0[62];
         rA2 = pA0[110];
         rA3 = pA0[158];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[15];
         rB0 = pB0[15];
         rA1 = pA0[63];
         rA2 = pA0[111];
         rA3 = pA0[159];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[16];
         rB0 = pB0[16];
         rA1 = pA0[64];
         rA2 = pA0[112];
         rA3 = pA0[160];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[17];
         rB0 = pB0[17];
         rA1 = pA0[65];
         rA2 = pA0[113];
         rA3 = pA0[161];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[18];
         rB0 = pB0[18];
         rA1 = pA0[66];
         rA2 = pA0[114];
         rA3 = pA0[162];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[19];
         rB0 = pB0[19];
         rA1 = pA0[67];
         rA2 = pA0[115];
         rA3 = pA0[163];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[20];
         rB0 = pB0[20];
         rA1 = pA0[68];
         rA2 = pA0[116];
         rA3 = pA0[164];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[21];
         rB0 = pB0[21];
         rA1 = pA0[69];
         rA2 = pA0[117];
         rA3 = pA0[165];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[22];
         rB0 = pB0[22];
         rA1 = pA0[70];
         rA2 = pA0[118];
         rA3 = pA0[166];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[23];
         rB0 = pB0[23];
         rA1 = pA0[71];
         rA2 = pA0[119];
         rA3 = pA0[167];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[24];
         rB0 = pB0[24];
         rA1 = pA0[72];
         rA2 = pA0[120];
         rA3 = pA0[168];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[25];
         rB0 = pB0[25];
         rA1 = pA0[73];
         rA2 = pA0[121];
         rA3 = pA0[169];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[26];
         rB0 = pB0[26];
         rA1 = pA0[74];
         rA2 = pA0[122];
         rA3 = pA0[170];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[27];
         rB0 = pB0[27];
         rA1 = pA0[75];
         rA2 = pA0[123];
         rA3 = pA0[171];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[28];
         rB0 = pB0[28];
         rA1 = pA0[76];
         rA2 = pA0[124];
         rA3 = pA0[172];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[29];
         rB0 = pB0[29];
         rA1 = pA0[77];
         rA2 = pA0[125];
         rA3 = pA0[173];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[30];
         rB0 = pB0[30];
         rA1 = pA0[78];
         rA2 = pA0[126];
         rA3 = pA0[174];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[31];
         rB0 = pB0[31];
         rA1 = pA0[79];
         rA2 = pA0[127];
         rA3 = pA0[175];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[32];
         rB0 = pB0[32];
         rA1 = pA0[80];
         rA2 = pA0[128];
         rA3 = pA0[176];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[33];
         rB0 = pB0[33];
         rA1 = pA0[81];
         rA2 = pA0[129];
         rA3 = pA0[177];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[34];
         rB0 = pB0[34];
         rA1 = pA0[82];
         rA2 = pA0[130];
         rA3 = pA0[178];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[35];
         rB0 = pB0[35];
         rA1 = pA0[83];
         rA2 = pA0[131];
         rA3 = pA0[179];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[36];
         rB0 = pB0[36];
         rA1 = pA0[84];
         rA2 = pA0[132];
         rA3 = pA0[180];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[37];
         rB0 = pB0[37];
         rA1 = pA0[85];
         rA2 = pA0[133];
         rA3 = pA0[181];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[38];
         rB0 = pB0[38];
         rA1 = pA0[86];
         rA2 = pA0[134];
         rA3 = pA0[182];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[39];
         rB0 = pB0[39];
         rA1 = pA0[87];
         rA2 = pA0[135];
         rA3 = pA0[183];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[40];
         rB0 = pB0[40];
         rA1 = pA0[88];
         rA2 = pA0[136];
         rA3 = pA0[184];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[41];
         rB0 = pB0[41];
         rA1 = pA0[89];
         rA2 = pA0[137];
         rA3 = pA0[185];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[42];
         rB0 = pB0[42];
         rA1 = pA0[90];
         rA2 = pA0[138];
         rA3 = pA0[186];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[43];
         rB0 = pB0[43];
         rA1 = pA0[91];
         rA2 = pA0[139];
         rA3 = pA0[187];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[44];
         rB0 = pB0[44];
         rA1 = pA0[92];
         rA2 = pA0[140];
         rA3 = pA0[188];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[45];
         rB0 = pB0[45];
         rA1 = pA0[93];
         rA2 = pA0[141];
         rA3 = pA0[189];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[46];
         rB0 = pB0[46];
         rA1 = pA0[94];
         rA2 = pA0[142];
         rA3 = pA0[190];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[47];
         rB0 = pB0[47];
         rA1 = pA0[95];
         rA2 = pA0[143];
         rA3 = pA0[191];
         rC0_0 += rA0 * rB0;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         pA0 += incAk;
         pB0 += incBk;
         *pC0 = rC0_0;
         pC0[1] = rC1_0;
         pC0[2] = rC2_0;
         pC0[3] = rC3_0;
         pC0 += incCm;
         pA0 += incAm;
         pB0 += incBm;
      }
      while(pA0 != stM);
      pC0 += incCn;
      pA0 += incAn;
      pB0 += incBn;
   }
   while(pB0 != stN);
}
#ifdef incAm
   #undef incAm
#endif
#ifdef incAn
   #undef incAn
#endif
#ifdef incAk
   #undef incAk
#endif
#ifdef incBm
   #undef incBm
#endif
#ifdef incBn
   #undef incBn
#endif
#ifdef incBk
   #undef incBk
#endif
#ifdef incCm
   #undef incCm
#endif
#ifdef incCn
   #undef incCn
#endif
#ifdef incCk
   #undef incCk
#endif
#ifdef Mb
   #undef Mb
#endif
#ifdef Nb
   #undef Nb
#endif
#ifdef Kb
   #undef Kb
#endif
void ATL_dJIK48x0x48TN48x48x0_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=48, NB=0, KB=48, 
 * lda=48, ldb=48, ldc=0, mu=4, nu=2, ku=48, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const int Nb = (N>>1)<<1;
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + 2304;
   const double *stN = B + (((Nb) << 5)+((Nb) << 4));
   #define incAk 48
   const int incAm = 144, incAn = -2304;
   #define incBk 48
   const int incBm = -48, incBn = 96;
   #define incCm 4
   const int incCn = (((ldc) << 1)) - 48;
   double *pC0=C, *pC1=pC0+(ldc);
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0, rA1, rA2, rA3;
   register double rB0, rB1;
   register double rC0_0, rC1_0, rC2_0, rC3_0, rC0_1, rC1_1, rC2_1, rC3_1;
   if (pB0 != stN)
   {
      do /* N-loop */
      {
         do /* M-loop */
         {
            rA0 = beta;
            rC0_0 = *pC0;
            rC0_0 *= rA0;
            rC1_0 = pC0[1];
            rC1_0 *= rA0;
            rC2_0 = pC0[2];
            rC2_0 *= rA0;
            rC3_0 = pC0[3];
            rC3_0 *= rA0;
            rC0_1 = *pC1;
            rC0_1 *= rA0;
            rC1_1 = pC1[1];
            rC1_1 *= rA0;
            rC2_1 = pC1[2];
            rC2_1 *= rA0;
            rC3_1 = pC1[3];
            rC3_1 *= rA0;
            rA0 = *pA0;
            rB0 = *pB0;
            rA1 = pA0[48];
            rA2 = pA0[96];
            rA3 = pA0[144];
            rB1 = pB0[48];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rA1 = pA0[49];
            rA2 = pA0[97];
            rA3 = pA0[145];
            rB1 = pB0[49];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rA1 = pA0[50];
            rA2 = pA0[98];
            rA3 = pA0[146];
            rB1 = pB0[50];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rA1 = pA0[51];
            rA2 = pA0[99];
            rA3 = pA0[147];
            rB1 = pB0[51];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rA1 = pA0[52];
            rA2 = pA0[100];
            rA3 = pA0[148];
            rB1 = pB0[52];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rA1 = pA0[53];
            rA2 = pA0[101];
            rA3 = pA0[149];
            rB1 = pB0[53];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rA1 = pA0[54];
            rA2 = pA0[102];
            rA3 = pA0[150];
            rB1 = pB0[54];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rA1 = pA0[55];
            rA2 = pA0[103];
            rA3 = pA0[151];
            rB1 = pB0[55];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[8];
            rB0 = pB0[8];
            rA1 = pA0[56];
            rA2 = pA0[104];
            rA3 = pA0[152];
            rB1 = pB0[56];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[9];
            rB0 = pB0[9];
            rA1 = pA0[57];
            rA2 = pA0[105];
            rA3 = pA0[153];
            rB1 = pB0[57];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[10];
            rB0 = pB0[10];
            rA1 = pA0[58];
            rA2 = pA0[106];
            rA3 = pA0[154];
            rB1 = pB0[58];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[11];
            rB0 = pB0[11];
            rA1 = pA0[59];
            rA2 = pA0[107];
            rA3 = pA0[155];
            rB1 = pB0[59];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[12];
            rB0 = pB0[12];
            rA1 = pA0[60];
            rA2 = pA0[108];
            rA3 = pA0[156];
            rB1 = pB0[60];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[13];
            rB0 = pB0[13];
            rA1 = pA0[61];
            rA2 = pA0[109];
            rA3 = pA0[157];
            rB1 = pB0[61];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[14];
            rB0 = pB0[14];
            rA1 = pA0[62];
            rA2 = pA0[110];
            rA3 = pA0[158];
            rB1 = pB0[62];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[15];
            rB0 = pB0[15];
            rA1 = pA0[63];
            rA2 = pA0[111];
            rA3 = pA0[159];
            rB1 = pB0[63];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[16];
            rB0 = pB0[16];
            rA1 = pA0[64];
            rA2 = pA0[112];
            rA3 = pA0[160];
            rB1 = pB0[64];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[17];
            rB0 = pB0[17];
            rA1 = pA0[65];
            rA2 = pA0[113];
            rA3 = pA0[161];
            rB1 = pB0[65];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[18];
            rB0 = pB0[18];
            rA1 = pA0[66];
            rA2 = pA0[114];
            rA3 = pA0[162];
            rB1 = pB0[66];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[19];
            rB0 = pB0[19];
            rA1 = pA0[67];
            rA2 = pA0[115];
            rA3 = pA0[163];
            rB1 = pB0[67];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[20];
            rB0 = pB0[20];
            rA1 = pA0[68];
            rA2 = pA0[116];
            rA3 = pA0[164];
            rB1 = pB0[68];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[21];
            rB0 = pB0[21];
            rA1 = pA0[69];
            rA2 = pA0[117];
            rA3 = pA0[165];
            rB1 = pB0[69];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[22];
            rB0 = pB0[22];
            rA1 = pA0[70];
            rA2 = pA0[118];
            rA3 = pA0[166];
            rB1 = pB0[70];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[23];
            rB0 = pB0[23];
            rA1 = pA0[71];
            rA2 = pA0[119];
            rA3 = pA0[167];
            rB1 = pB0[71];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[24];
            rB0 = pB0[24];
            rA1 = pA0[72];
            rA2 = pA0[120];
            rA3 = pA0[168];
            rB1 = pB0[72];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[25];
            rB0 = pB0[25];
            rA1 = pA0[73];
            rA2 = pA0[121];
            rA3 = pA0[169];
            rB1 = pB0[73];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[26];
            rB0 = pB0[26];
            rA1 = pA0[74];
            rA2 = pA0[122];
            rA3 = pA0[170];
            rB1 = pB0[74];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[27];
            rB0 = pB0[27];
            rA1 = pA0[75];
            rA2 = pA0[123];
            rA3 = pA0[171];
            rB1 = pB0[75];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[28];
            rB0 = pB0[28];
            rA1 = pA0[76];
            rA2 = pA0[124];
            rA3 = pA0[172];
            rB1 = pB0[76];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[29];
            rB0 = pB0[29];
            rA1 = pA0[77];
            rA2 = pA0[125];
            rA3 = pA0[173];
            rB1 = pB0[77];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[30];
            rB0 = pB0[30];
            rA1 = pA0[78];
            rA2 = pA0[126];
            rA3 = pA0[174];
            rB1 = pB0[78];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[31];
            rB0 = pB0[31];
            rA1 = pA0[79];
            rA2 = pA0[127];
            rA3 = pA0[175];
            rB1 = pB0[79];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[32];
            rB0 = pB0[32];
            rA1 = pA0[80];
            rA2 = pA0[128];
            rA3 = pA0[176];
            rB1 = pB0[80];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[33];
            rB0 = pB0[33];
            rA1 = pA0[81];
            rA2 = pA0[129];
            rA3 = pA0[177];
            rB1 = pB0[81];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[34];
            rB0 = pB0[34];
            rA1 = pA0[82];
            rA2 = pA0[130];
            rA3 = pA0[178];
            rB1 = pB0[82];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[35];
            rB0 = pB0[35];
            rA1 = pA0[83];
            rA2 = pA0[131];
            rA3 = pA0[179];
            rB1 = pB0[83];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[36];
            rB0 = pB0[36];
            rA1 = pA0[84];
            rA2 = pA0[132];
            rA3 = pA0[180];
            rB1 = pB0[84];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[37];
            rB0 = pB0[37];
            rA1 = pA0[85];
            rA2 = pA0[133];
            rA3 = pA0[181];
            rB1 = pB0[85];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[38];
            rB0 = pB0[38];
            rA1 = pA0[86];
            rA2 = pA0[134];
            rA3 = pA0[182];
            rB1 = pB0[86];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[39];
            rB0 = pB0[39];
            rA1 = pA0[87];
            rA2 = pA0[135];
            rA3 = pA0[183];
            rB1 = pB0[87];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[40];
            rB0 = pB0[40];
            rA1 = pA0[88];
            rA2 = pA0[136];
            rA3 = pA0[184];
            rB1 = pB0[88];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[41];
            rB0 = pB0[41];
            rA1 = pA0[89];
            rA2 = pA0[137];
            rA3 = pA0[185];
            rB1 = pB0[89];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[42];
            rB0 = pB0[42];
            rA1 = pA0[90];
            rA2 = pA0[138];
            rA3 = pA0[186];
            rB1 = pB0[90];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[43];
            rB0 = pB0[43];
            rA1 = pA0[91];
            rA2 = pA0[139];
            rA3 = pA0[187];
            rB1 = pB0[91];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[44];
            rB0 = pB0[44];
            rA1 = pA0[92];
            rA2 = pA0[140];
            rA3 = pA0[188];
            rB1 = pB0[92];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[45];
            rB0 = pB0[45];
            rA1 = pA0[93];
            rA2 = pA0[141];
            rA3 = pA0[189];
            rB1 = pB0[93];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[46];
            rB0 = pB0[46];
            rA1 = pA0[94];
            rA2 = pA0[142];
            rA3 = pA0[190];
            rB1 = pB0[94];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[47];
            rB0 = pB0[47];
            rA1 = pA0[95];
            rA2 = pA0[143];
            rA3 = pA0[191];
            rB1 = pB0[95];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            pA0 += incAk;
            pB0 += incBk;
            *pC0 = rC0_0;
            pC0[1] = rC1_0;
            pC0[2] = rC2_0;
            pC0[3] = rC3_0;
            *pC1 = rC0_1;
            pC1[1] = rC1_1;
            pC1[2] = rC2_1;
            pC1[3] = rC3_1;
            pC0 += incCm;
            pC1 += incCm;
            pA0 += incAm;
            pB0 += incBm;
         }
         while(pA0 != stM);
         pC0 += incCn;
         pC1 += incCn;
         pA0 += incAn;
         pB0 += incBn;
      }
      while(pB0 != stN);
   }
   if (k=N-Nb)
      ATL_dJIK48x0x48TN4x1x48_a1_bX(48, k, 48, alpha, ca, lda, cb + (((Nb) << 5)+((Nb) << 4)), ldb, beta, cc + (Nb*ldc), ldc);
}
#ifdef incAm
   #undef incAm
#endif
#ifdef incAn
   #undef incAn
#endif
#ifdef incAk
   #undef incAk
#endif
#ifdef incBm
   #undef incBm
#endif
#ifdef incBn
   #undef incBn
#endif
#ifdef incBk
   #undef incBk
#endif
#ifdef incCm
   #undef incCm
#endif
#ifdef incCn
   #undef incCn
#endif
#ifdef incCk
   #undef incCk
#endif
#ifdef Mb
   #undef Mb
#endif
#ifdef Nb
   #undef Nb
#endif
#ifdef Kb
   #undef Kb
#endif
