#include	<math.h>
#include	<stdio.h>

main()
{

   int i,N;
   double x,r,dr,phi1,phi2;
   double phi0,d0,d1,dc,ma,mc;
   double c0,c1,c2,c3;

   phi0 = 8.18555;
   ma = 3.30304;
   mc = 8.6655;
   d0 = 1.64;
   dc = 2.1052;
   d1 = 2.57;

   c0 =  2.2504290109e-8;
   c1 = -1.4408640561e-6;
   c2 =  2.1043303374e-5;
   c3 =  6.6024390226e-5;

   dr = 0.05;
   N = 1000;
   for (i=1; i<=N; ++i)
   {
      r = i*dr;
      phi1 = phi0*pow((d0/r),ma)* exp(-ma*pow(r/dc,mc) + ma*pow(d0/dc,mc));
      x = r-d1;
      phi2 = c0 + c1*x + c2*x*x + c3*x*x*x;
      printf("%le   %le %le\n",r,phi1,phi2);
   }
}
