/* PGrid.C
   Author - Eric Bylaska

	this class is used for defining 3d parallel maps
*/

#include	"Parallel.h"
#include	"PGrid.h"
#include	"control.h"
#include	"lattice.h"
#include	"blas.h"


/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/

PGrid::PGrid(Parallel *inparall) : d3db(inparall,control_mapping(),control_ngrid(0),control_ngrid(1),control_ngrid(2))
{
   int i,j,k,nxh,nyh,nzh,p,indx,k1,k2,k3,nb;
   int nwave_in[2],nwave_out[2];
   double *G1, *G2,*G3;
   double gx,gy,gz,gg,ggcut,eps;

   eps = 1.0e-12;
   Garray = new double [3*nfft3d];
   G1 = Garray; 
   G2 = &Garray[nfft3d]; 
   G3 = &Garray[2*nfft3d];
   nxh = nx/2;
   nyh = ny/2;
   nzh = nz/2;
   for (k3 = (-nzh+1); k3<= nzh; ++k3)
   for (k2 = (-nyh+1); k2<= nyh; ++k2)
   for (k1 = 0;        k1<= nxh; ++k1)
   {
      gx = k1*lattice_unitg(0,0) + k2*lattice_unitg(0,1) + k3*lattice_unitg(0,2);
      gy = k1*lattice_unitg(1,0) + k2*lattice_unitg(1,1) + k3*lattice_unitg(1,2);
      gz = k1*lattice_unitg(2,0) + k2*lattice_unitg(2,1) + k3*lattice_unitg(2,2);
      i=k1; if (i < 0) i = i + nx;
      j=k2; if (j < 0) j = j + ny;
      k=k3; if (k < 0) k = k + nz;

      indx = ijktoindex(i,j,k);
      p    = ijktop(i,j,k);
      if (p==parall->taskid_i())
      {
         G1[indx] = gx;
         G2[indx] = gy;
         G3[indx] = gz;
      }

   }
   masker[0] = new int [2*nfft3d];
   masker[1] = &(masker[0][nfft3d]);
   for (int k=0; k<(2*nfft3d); ++k) 
      masker[0][k] = 1;

   for (nb=0; nb<=1; ++nb)
   {
      nwave[nb] = 0;
      if (nb==0) 
         ggcut = lattice_eggcut();
      else
         ggcut = lattice_wggcut();

      for (k3 = (-nzh+1); k3<  nzh; ++k3)
      for (k2 = (-nyh+1); k2<  nyh; ++k2)
      for (k1 = 0;        k1<  nxh; ++k1)
      {
         gx = k1*lattice_unitg(0,0) + k2*lattice_unitg(0,1) + k3*lattice_unitg(0,2);
         gy = k1*lattice_unitg(1,0) + k2*lattice_unitg(1,1) + k3*lattice_unitg(1,2);
         gz = k1*lattice_unitg(2,0) + k2*lattice_unitg(2,1) + k3*lattice_unitg(2,2);
         i=k1; if (i < 0) i = i + nx;
         j=k2; if (j < 0) j = j + ny;
         k=k3; if (k < 0) k = k + nz;

         indx = ijktoindex(i,j,k);
         p    = ijktop(i,j,k);
         if (p==parall->taskid_i())
         {
            gg = gx*gx + gy*gy + gz*gz;
            gg = gg-ggcut;
            if (gg < (-eps))
            {
               masker[nb][indx] = 0;
               ++nwave[nb] ;
            }
         }
      }
      nwave_entire[nb] = nwave[nb];
      nwave_entire[nb] = parall->ISumAll(1,nwave_entire[nb]);
   }


   packarray[0] = new int [2*nfft3d];
   packarray[1] = &(packarray[0][nfft3d]);

   for (nb=0; nb<=1; ++nb)
   {
      nida[nb]  = 0;
      nidb2[nb] = 0;

      /* k=(0,0,0)  */
      k1=0;
      k2=0;
      k3=0;
      indx = ijktoindex(k1,k2,k3);
      p    = ijktop(k1,k2,k3);
      if (p==parall->taskid_i())
         if (!masker[nb][indx])
         {
            packarray[nb][nida[nb]] = indx;
            ++nida[nb];
         }

      /* k=(0,0,k3) - neglect (0,0,-k3) points */
      for (k=1; k<=(nzh-1); ++k)
      {
         k1=0;
         k2=0;
         k3=k;
         indx = ijktoindex(k1,k2,k3);
         p    = ijktop(k1,k2,k3);
         if (p==parall->taskid_i())
            if (!masker[nb][indx])
            {
               packarray[nb][nida[nb]+nidb2[nb]] = indx;
               ++nidb2[nb];
            }
      }

      /* k=(0,k2,k3) - neglect (0,-k2, -k3) points */
      for (k=(-nzh+1); k<=(nzh-1); ++k)
      for (j=1;        j<=(nyh-1); ++j)
      {
         k1=0;
         k2=j;
         k3=k;
         if (k3 < 0) k3 = k3 + nz;
         indx = ijktoindex(k1,k2,k3);
         p    = ijktop(k1,k2,k3);
         if (p==parall->taskid_i())
            if (!masker[nb][indx])
            {
               packarray[nb][nida[nb]+nidb2[nb]] = indx;
               ++nidb2[nb];
            }
      }


      /* k=(k1,k2,k3) */
      for (k=(-nzh+1); k<=(nzh-1); ++k)
      for (j=(-nyh+1); j<=(nyh-1); ++j)
      for (i=1;        i<=(nxh-1); ++i)
      {
         k1=i;
         k2=j;
         k3=k;
         if (k2 < 0) k2 = k2 + ny;
         if (k3 < 0) k3 = k3 + nz;
         indx = ijktoindex(k1,k2,k3);
         p    = ijktop(k1,k2,k3);
         if (p==parall->taskid_i())
            if (!masker[nb][indx])
            {
               packarray[nb][nida[nb]+nidb2[nb]] = indx;
               ++nidb2[nb];
            }
      }

   }

   nwave_in[0] = nida[0] + nidb2[0];
   nwave_in[1] = nida[1] + nidb2[1];


   if (control_balance()) 
   {
      balanced = 1;
      mybalance = new Balance(parall,2,nwave_in,nwave_out);
   }
   else
   {
      balanced = 0;
      nwave_out[0] = nwave_in[0];
      nwave_out[1] = nwave_in[1];
   }
   nidb[0] = nidb2[0] + (nwave_out[0] - nwave_in[0]);
   nidb[1] = nidb2[1] + (nwave_out[1] - nwave_in[1]);

   nwave_all[0] = nida[0] + nidb2[0];
   nwave_all[1] = nida[1] + nidb2[1];
   parall->Vector_ISumAll(1,2,nwave_all);
}


void c_indexcopy(const int n, const int *indx, double *A, double *B)
{
   int ii,jj;
   ii = 0;
   for (int i=0; i<n; ++i)
   {
      jj      = 2*indx[i];
      B[ii]   = A[jj];
      B[ii+1] = A[jj+1];
      ii +=2;
   }
}


void t_indexcopy(const int n, const int *indx, double *A, double *B)
{
   for (int i=0; i<n; ++i)
      B[i] = A[indx[i]];
}

void PGrid::c_pack(const int nb, double *a)
{
   int one      = 1;
   int zero     = 0;
   double rzero = 0.0;
   double *tmp = new double [2*nfft3d];

   dcopy_(&n2ft3d,a,&one,tmp,&one);
   dcopy_(&n2ft3d,&rzero,&zero,a,&one);
   c_indexcopy(nida[nb]+nidb2[nb],packarray[nb],tmp,a);

   //if (balanced)
   //   mybalance->c_balance()

   delete [] tmp;
}

void PGrid::t_pack(const int nb, double *a)
{
   int one      = 1;
   int zero     = 0;
   double rzero = 0.0;
   double *tmp = new double [nfft3d];

   dcopy_(&nfft3d,a,&one,tmp,&one);
   dcopy_(&nfft3d,&rzero,&zero,a,&one);
   t_indexcopy(nida[nb]+nidb2[nb],packarray[nb],tmp,a);

   //if (balanced)
   //   mybalance->t_balance()

   delete [] tmp;
}

void PGrid::tt_pack_copy(const int nb, double *a, double *b)
{
   int one = 1;
   int ng  = nida[nb]+nidb[nb];
   dcopy_(&ng,a,&one,b,&one);
}


