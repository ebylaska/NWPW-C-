

//#include	"float.h"

/**************************************
 *                                    *
 *            eigsrt                  *
 *                                    *
 **************************************/

void eigsrt(double *D, double *V, int n)
{
   int i,j,k;
   double p;

   for (i=0; i<(n-1); ++i)
   {
      k = i;
      p = D[i];
      for(j=i+1; j<n; ++j)
         if (D[j]>=p)
         {
            k = j;
            p = D[j];
         }

      if (k!=i)
      {
         D[k] = D[i];
         D[i] = p;
         for (j=0; j<n; ++j)
         {
            p = V[j+i*n];
            V[j+i*n] = V[j+k*n];
            V[j+k*n] = p;
         }
      }
   }
}

