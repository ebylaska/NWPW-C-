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
static void ATL_zJIK52x0x52TN4x1x52_a1_b1
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=52, NB=0, KB=52, 
 * lda=52, ldb=52, ldc=0, mu=4, nu=1, ku=52, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   #define Nb N
   const double *stM = A + 2704;
   const double *stN = B + (52*(Nb));
   #define incAk 52
   const int incAm = 156, incAn = -2704;
   #define incBk 52
   const int incBm = -52, incBn = 52;
   #define incCm 8
   const int incCn = (((ldc) << 1)) - 104;
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
         rC0_0 = *pC0;
         rC1_0 = pC0[2];
         rC2_0 = pC0[4];
         rC3_0 = pC0[6];
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[52];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[104];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[156];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[1];
         rB0 = pB0[1];
         rA1 = pA0[53];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[105];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[157];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[2];
         rB0 = pB0[2];
         rA1 = pA0[54];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[106];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[158];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[3];
         rB0 = pB0[3];
         rA1 = pA0[55];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[107];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[159];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[4];
         rB0 = pB0[4];
         rA1 = pA0[56];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[108];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[160];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[5];
         rB0 = pB0[5];
         rA1 = pA0[57];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[109];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[161];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[6];
         rB0 = pB0[6];
         rA1 = pA0[58];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[110];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[162];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[7];
         rB0 = pB0[7];
         rA1 = pA0[59];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[111];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[163];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[8];
         rB0 = pB0[8];
         rA1 = pA0[60];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[112];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[164];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[9];
         rB0 = pB0[9];
         rA1 = pA0[61];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[113];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[165];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[10];
         rB0 = pB0[10];
         rA1 = pA0[62];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[114];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[166];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[11];
         rB0 = pB0[11];
         rA1 = pA0[63];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[115];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[167];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[12];
         rB0 = pB0[12];
         rA1 = pA0[64];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[116];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[168];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[13];
         rB0 = pB0[13];
         rA1 = pA0[65];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[117];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[169];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[14];
         rB0 = pB0[14];
         rA1 = pA0[66];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[118];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[170];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[15];
         rB0 = pB0[15];
         rA1 = pA0[67];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[119];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[171];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[16];
         rB0 = pB0[16];
         rA1 = pA0[68];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[120];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[172];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[17];
         rB0 = pB0[17];
         rA1 = pA0[69];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[121];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[173];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[18];
         rB0 = pB0[18];
         rA1 = pA0[70];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[122];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[174];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[19];
         rB0 = pB0[19];
         rA1 = pA0[71];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[123];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[175];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[20];
         rB0 = pB0[20];
         rA1 = pA0[72];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[124];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[176];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[21];
         rB0 = pB0[21];
         rA1 = pA0[73];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[125];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[177];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[22];
         rB0 = pB0[22];
         rA1 = pA0[74];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[126];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[178];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[23];
         rB0 = pB0[23];
         rA1 = pA0[75];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[127];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[179];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[24];
         rB0 = pB0[24];
         rA1 = pA0[76];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[128];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[180];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[25];
         rB0 = pB0[25];
         rA1 = pA0[77];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[129];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[181];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[26];
         rB0 = pB0[26];
         rA1 = pA0[78];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[130];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[182];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[27];
         rB0 = pB0[27];
         rA1 = pA0[79];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[131];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[183];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[28];
         rB0 = pB0[28];
         rA1 = pA0[80];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[132];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[184];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[29];
         rB0 = pB0[29];
         rA1 = pA0[81];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[133];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[185];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[30];
         rB0 = pB0[30];
         rA1 = pA0[82];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[134];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[186];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[31];
         rB0 = pB0[31];
         rA1 = pA0[83];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[135];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[187];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[32];
         rB0 = pB0[32];
         rA1 = pA0[84];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[136];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[188];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[33];
         rB0 = pB0[33];
         rA1 = pA0[85];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[137];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[189];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[34];
         rB0 = pB0[34];
         rA1 = pA0[86];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[138];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[190];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[35];
         rB0 = pB0[35];
         rA1 = pA0[87];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[139];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[191];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[36];
         rB0 = pB0[36];
         rA1 = pA0[88];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[140];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[192];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[37];
         rB0 = pB0[37];
         rA1 = pA0[89];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[141];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[193];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[38];
         rB0 = pB0[38];
         rA1 = pA0[90];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[142];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[194];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[39];
         rB0 = pB0[39];
         rA1 = pA0[91];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[143];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[195];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[40];
         rB0 = pB0[40];
         rA1 = pA0[92];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[144];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[196];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[41];
         rB0 = pB0[41];
         rA1 = pA0[93];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[145];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[197];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[42];
         rB0 = pB0[42];
         rA1 = pA0[94];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[146];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[198];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[43];
         rB0 = pB0[43];
         rA1 = pA0[95];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[147];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[199];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[44];
         rB0 = pB0[44];
         rA1 = pA0[96];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[148];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[200];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[45];
         rB0 = pB0[45];
         rA1 = pA0[97];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[149];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[201];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[46];
         rB0 = pB0[46];
         rA1 = pA0[98];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[150];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[202];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[47];
         rB0 = pB0[47];
         rA1 = pA0[99];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[151];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[203];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[48];
         rB0 = pB0[48];
         rA1 = pA0[100];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[152];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[204];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[49];
         rB0 = pB0[49];
         rA1 = pA0[101];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[153];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[205];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[50];
         rB0 = pB0[50];
         rA1 = pA0[102];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[154];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[206];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rA0 = pA0[51];
         rB0 = pB0[51];
         rA1 = pA0[103];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[155];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[207];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         pA0 += incAk;
         pB0 += incBk;
         *pC0 = rC0_0;
         pC0[2] = rC1_0;
         pC0[4] = rC2_0;
         pC0[6] = rC3_0;
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
void ATL_zJIK52x0x52TN52x52x0_a1_b1
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=52, NB=0, KB=52, 
 * lda=52, ldb=52, ldc=0, mu=4, nu=2, ku=52, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const int Nb = (N>>1)<<1;
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + 2704;
   const double *stN = B + (52*(Nb));
   #define incAk 52
   const int incAm = 156, incAn = -2704;
   #define incBk 52
   const int incBm = -52, incBn = 104;
   #define incCm 8
   const int incCn = (((ldc) << 2)) - 104;
   double *pC0=C, *pC1=pC0+(((ldc) << 1));
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
            rC0_0 = *pC0;
            rC1_0 = pC0[2];
            rC2_0 = pC0[4];
            rC3_0 = pC0[6];
            rC0_1 = *pC1;
            rC1_1 = pC1[2];
            rC2_1 = pC1[4];
            rC3_1 = pC1[6];
            rA0 = *pA0;
            rB0 = *pB0;
            rA1 = pA0[52];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[104];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[156];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[52];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rA1 = pA0[53];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[105];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[157];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[53];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rA1 = pA0[54];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[106];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[158];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[54];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rA1 = pA0[55];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[107];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[159];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[55];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rA1 = pA0[56];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[108];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[160];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[56];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rA1 = pA0[57];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[109];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[161];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[57];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rA1 = pA0[58];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[110];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[162];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[58];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rA1 = pA0[59];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[111];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[163];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[59];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[8];
            rB0 = pB0[8];
            rA1 = pA0[60];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[112];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[164];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[60];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[9];
            rB0 = pB0[9];
            rA1 = pA0[61];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[113];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[165];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[61];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[10];
            rB0 = pB0[10];
            rA1 = pA0[62];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[114];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[166];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[62];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[11];
            rB0 = pB0[11];
            rA1 = pA0[63];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[115];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[167];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[63];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[12];
            rB0 = pB0[12];
            rA1 = pA0[64];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[116];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[168];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[64];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[13];
            rB0 = pB0[13];
            rA1 = pA0[65];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[117];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[169];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[65];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[14];
            rB0 = pB0[14];
            rA1 = pA0[66];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[118];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[170];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[66];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[15];
            rB0 = pB0[15];
            rA1 = pA0[67];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[119];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[171];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[67];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[16];
            rB0 = pB0[16];
            rA1 = pA0[68];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[120];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[172];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[68];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[17];
            rB0 = pB0[17];
            rA1 = pA0[69];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[121];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[173];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[69];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[18];
            rB0 = pB0[18];
            rA1 = pA0[70];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[122];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[174];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[70];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[19];
            rB0 = pB0[19];
            rA1 = pA0[71];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[123];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[175];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[71];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[20];
            rB0 = pB0[20];
            rA1 = pA0[72];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[124];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[176];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[72];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[21];
            rB0 = pB0[21];
            rA1 = pA0[73];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[125];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[177];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[73];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[22];
            rB0 = pB0[22];
            rA1 = pA0[74];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[126];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[178];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[74];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[23];
            rB0 = pB0[23];
            rA1 = pA0[75];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[127];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[179];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[75];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[24];
            rB0 = pB0[24];
            rA1 = pA0[76];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[128];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[180];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[76];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[25];
            rB0 = pB0[25];
            rA1 = pA0[77];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[129];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[181];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[77];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[26];
            rB0 = pB0[26];
            rA1 = pA0[78];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[130];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[182];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[78];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[27];
            rB0 = pB0[27];
            rA1 = pA0[79];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[131];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[183];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[79];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[28];
            rB0 = pB0[28];
            rA1 = pA0[80];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[132];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[184];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[80];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[29];
            rB0 = pB0[29];
            rA1 = pA0[81];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[133];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[185];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[81];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[30];
            rB0 = pB0[30];
            rA1 = pA0[82];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[134];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[186];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[82];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[31];
            rB0 = pB0[31];
            rA1 = pA0[83];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[135];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[187];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[83];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[32];
            rB0 = pB0[32];
            rA1 = pA0[84];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[136];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[188];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[84];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[33];
            rB0 = pB0[33];
            rA1 = pA0[85];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[137];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[189];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[85];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[34];
            rB0 = pB0[34];
            rA1 = pA0[86];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[138];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[190];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[86];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[35];
            rB0 = pB0[35];
            rA1 = pA0[87];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[139];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[191];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[87];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[36];
            rB0 = pB0[36];
            rA1 = pA0[88];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[140];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[192];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[88];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[37];
            rB0 = pB0[37];
            rA1 = pA0[89];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[141];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[193];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[89];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[38];
            rB0 = pB0[38];
            rA1 = pA0[90];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[142];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[194];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[90];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[39];
            rB0 = pB0[39];
            rA1 = pA0[91];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[143];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[195];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[91];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[40];
            rB0 = pB0[40];
            rA1 = pA0[92];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[144];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[196];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[92];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[41];
            rB0 = pB0[41];
            rA1 = pA0[93];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[145];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[197];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[93];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[42];
            rB0 = pB0[42];
            rA1 = pA0[94];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[146];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[198];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[94];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[43];
            rB0 = pB0[43];
            rA1 = pA0[95];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[147];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[199];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[95];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[44];
            rB0 = pB0[44];
            rA1 = pA0[96];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[148];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[200];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[96];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[45];
            rB0 = pB0[45];
            rA1 = pA0[97];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[149];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[201];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[97];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[46];
            rB0 = pB0[46];
            rA1 = pA0[98];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[150];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[202];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[98];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[47];
            rB0 = pB0[47];
            rA1 = pA0[99];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[151];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[203];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[99];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[48];
            rB0 = pB0[48];
            rA1 = pA0[100];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[152];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[204];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[100];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[49];
            rB0 = pB0[49];
            rA1 = pA0[101];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[153];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[205];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[101];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[50];
            rB0 = pB0[50];
            rA1 = pA0[102];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[154];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[206];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[102];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            rA0 = pA0[51];
            rB0 = pB0[51];
            rA1 = pA0[103];
            rC0_0 += rA0 * rB0;
            rA2 = pA0[155];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[207];
            rC2_0 += rA2 * rB0;
            rB1 = pB0[103];
            rC3_0 += rA3 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC2_1 += rA2 * rB1;
            rC3_1 += rA3 * rB1;
            pA0 += incAk;
            pB0 += incBk;
            *pC0 = rC0_0;
            pC0[2] = rC1_0;
            pC0[4] = rC2_0;
            pC0[6] = rC3_0;
            *pC1 = rC0_1;
            pC1[2] = rC1_1;
            pC1[4] = rC2_1;
            pC1[6] = rC3_1;
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
      ATL_zJIK52x0x52TN4x1x52_a1_b1(52, k, 52, alpha, ca, lda, cb + (52*(Nb)), ldb, beta, cc + (((Nb*ldc) << 1)), ldc);
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
#ifdef ATL_UCLEANN
#define ATL_zpNBmm_b1 ATL_zgpNBmm_b1
#endif

void ATL_zJIK52x0x52TN52x52x0_a1_bX(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);
void ATL_zJIK52x0x52TN52x52x0_a1_bX(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);
void ATL_zJIK52x0x52TN52x52x0_a1_b1(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);

void ATL_zpNBmm_b1(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc)
{
   ATL_zJIK52x0x52TN52x52x0_a1_bX(M, N, K, alpha, A, lda, B, ldb, -beta, C, ldc);
   ATL_zJIK52x0x52TN52x52x0_a1_bX(M, N, K, alpha, A, lda, B+N*ldb, ldb, beta, C+1, ldc);
   ATL_zJIK52x0x52TN52x52x0_a1_bX(M, N, K, alpha, A+M*lda, lda, B+N*ldb, ldb, -1.0, C, ldc);
   ATL_zJIK52x0x52TN52x52x0_a1_b1(M, N, K, alpha, A+M*lda, lda, B, ldb, 1.0, C+1, ldc);
}
