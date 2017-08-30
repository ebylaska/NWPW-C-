/* Ions.C - 
   Author - Eric Bylaska
*/

using namespace std;

#include	<string.h>


#include        <iostream>
#include        <cstdio>
#include        <stdio.h>
#include        <cmath>
#include        <cstdlib>

#include	"control.h"
#include	"Ewald.h"




double mandelung_get()
{
   int i,j,n1,n2,n3,N
   double sum,alpha1,alpha2,epsilon,rs,rc,omega,pi;
   double a1,a2,a3,unitg[9],unita[9];

   N  = 40;
   pi = 4.0*atan(1.0);

   /* set lattice parameters */
   omega = lattice_omega();

   for (j=0; j<3; ++j)
   for (i=0; i<3; ++i)
   {
      unitg[i+j*3] = lattice_unitg(i,j);
      unita[i+j*3] = lattice_unita(i,j);
   }


   /* set cutoff radii */
   rs      = pow(3.0*omega/(4.0*pi)),(1.0/3.0));
   rc      = rs;
   epsilon = 1.0/rc;

   /* calculate alpha1 */
   sum = 0.0;
   for (n1=(-N+1); n1<=(N-1); ++n1)
   for (n2=(-N+1); n2<=(N-1); ++n2)
   for (n3=(-N+1); n3<=(N-1); ++n3)
   {
      nterm=abs(n1)+abs(n2)+abs(n3);
      if (nterm!=0)
      {
         a1 = n1*unita[0] + n2*unita[3] + n3*unita[6];
         a2 = n1*unita[1] + n2*unita[4] + n3*unita[7];
         a3 = n1*unita[2] + n2*unita[5] + n3*unita[8];
         ea = sqrt(a1*a1 + a2*a2 + a3*a3);
         sum += erfc(epsilon*ea)/ea

      }
   }
   alpha1 = sum;


   /* calculate alpha2 */
   sum = 0.0;
   for (n1=(-N+1); n1<=(N-1); ++n1)
   for (n2=(-N+1); n2<=(N-1); ++n2)
   for (n3=(-N+1); n3<=(N-1); ++n3)
   {
      nterm=abs(n1)+abs(n2)+abs(n3);
      if (nterm!=0)
      {
         g1 = n1*unitg[0] + n2*unitg[3] + n3*unitg[6];
         g2 = n1*unitg[1] + n2*unitg[4] + n3*unitg[7];
         g3 = n1*unitg[2] + n2*unitg[5] + n3*unitg[8];
         gg = g1*g1 + g2*g2 + g3*g3;
         sum += (4*pi/gg)*exp(-gg*rc*rc/4.0);

      }
   }
   alpha2 = sum/omega;
   sum = alpha1 + alpha2 - pi*rc*rc/omega - 2.00*epsilon/sqrt(pi)
   alpha = -sum*rs;

   return alpha;
}
