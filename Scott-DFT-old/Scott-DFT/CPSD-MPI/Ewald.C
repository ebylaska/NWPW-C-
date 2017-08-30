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


/* Constructors */

/*********************************
 *                               *
 *          Ewald::Ewald         *
 *                               *
 *********************************/
Ewald::Ewald(Parallel *inparall, Ion *inion)
{
   int i,j,k,k1,k2,k3;
   int enxh,enyh,enzh;
   int tnp,tid,dutask;
   double g1,g2,g3,gg1,gg2,gg3,gg,ggcut;
   double pi,pi4,rs,w;
   double eps=1.0e-12;

   ewaldparall = inparall;
   ewaldion    = inion;
   tnp = ewaldparall->np();
   tid = ewaldparall->taskid();

   for (j=0; j<3; ++j)
   for (i=0; i<3; ++i)
   {
      unitg[i+j*3] = lattice_unitg(i,j);
      unita[i+j*3] = lattice_unita(i,j);
   }
      
   enx=control_ewald_ngrid(0);
   eny=control_ewald_ngrid(1);
   enz=control_ewald_ngrid(2);
   enxh=enx/2;
   enyh=eny/2;
   enzh=enz/2;


   /* determine ggcut */
   g1 = unitg[0]*(enxh);
   g2 = unitg[1]*(enxh);
   g3 = unitg[2]*(enxh);
   gg1 = g1*g1 + g2*g2 + g3*g3;

   g1 = unitg[3]*(enyh);
   g2 = unitg[4]*(enyh);
   g3 = unitg[5]*(enyh);
   gg2 = g1*g1 + g2*g2 + g3*g3;

   g1 = unitg[6]*(enzh);
   g2 = unitg[7]*(enzh);
   g3 = unitg[8]*(enzh);
   gg3 = g1*g1 + g2*g2 + g3*g3;

   ggcut = gg1;
   if (gg2<ggcut) ggcut = gg2;
   if (gg3<ggcut) ggcut = gg3;
   if ((2.0*control_ecut())<ggcut) 
      ggcut=2.0*control_ecut();
   eecut = 0.5*ggcut;


   /* determine enpack */
   dutask= 0;
   enpack= 0;
   enida = 0;
   k1 = 0;
   k2 = 0;
   k3 = 0;
   g1=k1*unitg[0]+k2*unitg[3]+k3*unitg[6];
   g2=k1*unitg[1]+k2*unitg[4]+k3*unitg[7];
   g3=k1*unitg[2]+k2*unitg[5]+k3*unitg[8];
   gg=g1*g1+g2*g2+g3*g3;
   if ((gg-ggcut)<(-eps))
   {
      if (dutask==tid)
      { 
         enpack++;
         enida++;
      }
      dutask = (dutask+1)%tnp;
   }
   for (k=1; k<enzh; ++k)
   {
      k1 = 0;
      k2 = 0;
      k3 = k;
      g1=k1*unitg[0]+k2*unitg[3]+k3*unitg[6];
      g2=k1*unitg[1]+k2*unitg[4]+k3*unitg[7];
      g3=k1*unitg[2]+k2*unitg[5]+k3*unitg[8];
      gg=g1*g1+g2*g2+g3*g3;
      if ((gg-ggcut)<(-eps))
      {
         if (dutask==tid) enpack++;
         dutask =(dutask+1)%tnp;
      }
   }

   for (k=(-enzh+1); k<enzh; ++k)
   for (j=1; j<enyh; ++j)
   {
      k1 = 0;
      k2 = j;
      k3 = k;
      g1=k1*unitg[0]+k2*unitg[3]+k3*unitg[6];
      g2=k1*unitg[1]+k2*unitg[4]+k3*unitg[7];
      g3=k1*unitg[2]+k2*unitg[5]+k3*unitg[8];
      gg=g1*g1+g2*g2+g3*g3;
      if ((gg-ggcut)<(-eps))
      {
         if (dutask==tid) enpack++;
         dutask = (dutask+1)%tnp;
      }
   }
   

   for (k=(-enzh+1); k<enzh; ++k)
   for (j=(-enyh+1); j<enyh; ++j)
   for (i=1; i<enxh; ++i)
   {
      k1=i;
      k2=j;
      k3=k;
      g1=k1*unitg[0]+k2*unitg[3]+k3*unitg[6];
      g2=k1*unitg[1]+k2*unitg[4]+k3*unitg[7];
      g3=k1*unitg[2]+k2*unitg[5]+k3*unitg[8];
      gg=g1*g1+g2*g2+g3*g3;

      if ((gg-ggcut)< (-eps))
      {
         if (dutask==tid) enpack++;
         dutask = (dutask+1)%tnp;
      }
   }
   enpack_all = ewaldparall->ISumAll(0,enpack);



   encut = control_ewald_ncut();
   enshl3d = (2*encut+1)*(2*encut+1)*(2*encut+1);
   ercut = control_ewald_rcut();
   pi  = 4.00*atan(1.0);
   pi4 = 4.00*pi;
   if (encut<=0) encut=1;
   if (ercut<=0.00)
   {
      rs = unita[0]*unita[0] + unita[1]*unita[1] + unita[2]*unita[2];
      rs = sqrt(rs);
      ercut=rs/pi;

      rs = unita[3]*unita[3] + unita[4]*unita[4] + unita[5]*unita[5];
      rs = sqrt(rs);
      w=rs/pi;
      if (w<ercut) ercut = w;

      rs = unita[6]*unita[6] + unita[7]*unita[7] + unita[8]*unita[8];
      rs = sqrt(rs);
      w=rs/pi;
      if (w<ercut) ercut = w;
   }
   w = 0.25*ercut*ercut;


   /* allocate memory */
   eG  = new double [3*enpack];
   vg  = new double [enpack];
   vcx  = new double [enpack];
   rcell = new double [enshl3d];
   ewx1 = new double [2*(ewaldion->nion)*enx];
   ewy1 = new double [2*(ewaldion->nion)*eny];
   ewz1 = new double [2*(ewaldion->nion)*enz];
   zv   = new double [ewaldion->nkatm];
   i_indx = new int[enpack];
   j_indx = new int[enpack];
   k_indx = new int[enpack];

}

double mandelung_get()
{
   double alpha;

      pi = 4.0d0*datan(1.0d0)

*     ***** set lattice parameters *****
      omega = lattice_omega()

   for (j=0; j<3; ++j)
   for (i=0; i<3; ++i)
   {
      unitg[i+j*3] = lattice_unitg(i,j);
      unita[i+j*3] = lattice_unita(i,j);
   }


*     ***** set cutoff radii ****
      rs      = (3.0d0*omega/(4.0d0*pi))**(1.0d0/3.0d0)
      rc      = rs
      epsilon = 1.0d0/rc

*     **** calculate alpha1 *****
      sum = 0.0d0
      do n1=(-N+1),(N-1)
      do n2=(-N+1),(N-1)
      do n3=(-N+1),(N-1)
         nterm=iabs(n1)+iabs(n2)+iabs(n3)
         if (nterm.ne.0) then
            a1 = n1*unita(1,1)
     >         + n2*unita(1,2)
     >         + n3*unita(1,3)

            a2 = n1*unita(2,1)
     >         + n2*unita(2,2)
     >         + n3*unita(2,3)

            a3 = n1*unita(3,1)
     >         + n2*unita(3,2)
     >         + n3*unita(3,3)

            ea = dsqrt(a1*a1 + a2*a2 + a3*a3)

            sum = sum + util_erfc(epsilon*ea)/ea

         end if
      end do
      end do
      end do
      alpha1 = sum

                
*     **** calculate alpha2 *****
      sum = 0.0d0
      do n1=(-N+1),(N-1)
      do n2=(-N+1),(N-1)
      do n3=(-N+1),(N-1) 
         nterm=iabs(n1)+iabs(n2)+iabs(n3)
         if (nterm.ne.0) then
            g1 = n1*unitg(1,1)
     >         + n2*unitg(1,2)
     >         + n3*unitg(1,3)

            g2 = n1*unitg(2,1)
     >         + n2*unitg(2,2)
     >         + n3*unitg(2,3)

            g3 = n1*unitg(3,1)
     >         + n2*unitg(3,2)
     >         + n3*unitg(3,3)

            gg  = g1*g1 + g2*g2 + g3*g3 
            sum = sum +  (4.0d0*pi/gg)* exp(-gg*rc*rc/4.0d0)
             
         end if
      end do
      end do
      end do 
      alpha2 = sum/omega



      sum = alpha1 + alpha2  
     >    - pi*rc*rc/omega - 2.0d0*epsilon/dsqrt(pi)

      alpha = -sum*rs



   return alpha;
}
