/* Pseudopotential.C - 
   Author - Eric Bylaska
*/

using namespace std;

#include	<string.h>
#include        <iostream>
#include        <cstdio>
#include        <stdio.h>
#include        <cmath>
#include        <cstdlib>

#include	"compressed_io.h"
#include	"Pseudopotential.h"


static void psp_read(PGrid& mygrid,
                     char *fname, 
                     char *comment,
                     int *psp_type,
                     int *version,
                     int *nfft,
                     double *unita,
                     char   *atom,
                     double *amass,
                     double *zv,
                     int *lmmax,
                     int *lmax,
                     int *locp,
                     int *nmax,
                     double **rc,
                     int *nprj,
                     int **n_projector,
                     int **l_projector,
                     int **m_projector,
                     double **Gijl,
                     int *semicore,
                     double *rcore,
                     double **ncore,
                     double *vl,
                     double **vnl)
{
   int i,nn;
   double   *tmp2,*prj;
   Parallel *parall = mygrid.parall;

   if (parall->is_master())
   {
      openfile(5,fname,"r");
      cread(5,comment,80);
      comment[79] = '\0';
      i = 78;
      while (comment[i] == ' ')
        comment[i--] = '\0';
 
      iread(5,psp_type,1);
      iread(5,version,1);
      iread(5,nfft,3);
      dread(5,unita,9);
      cread(5,atom,2);
      dread(5,amass,1);
      dread(5,zv,1);
      iread(5,lmax,1);
      iread(5,locp,1);
      iread(5,nmax,1);
   }
   parall->Brdcst_cValues(0,0,80,comment);
   parall->Brdcst_iValue(0,0,psp_type);
   parall->Brdcst_iValue(0,0,version);
   parall->Brdcst_iValues(0,0,3,nfft);
   parall->Brdcst_Values(0,0,9,unita);
   parall->Brdcst_cValues(0,0,2,atom);
   parall->Brdcst_Values(0,0,1,amass);
   parall->Brdcst_Values(0,0,1,zv);
   parall->Brdcst_iValue(0,0,lmax);
   parall->Brdcst_iValue(0,0,locp);
   parall->Brdcst_iValue(0,0,nmax);
   *lmmax=((*lmax)+1)*((*lmax)+1) - (2*(*locp)+1);

   *rc = new double[*lmax+1];
   if (parall->is_master())
   {
      dread(5,*rc,*lmax+1);
      iread(5,nprj,1);
   }
   parall->Brdcst_Values(0,0,*lmax+1,*rc);
   parall->Brdcst_iValue(0,0,nprj);
   if (*nprj > 0) 
   {
      *n_projector = new int[*nprj];
      *l_projector = new int[*nprj];
      *m_projector = new int[*nprj];
      if (parall->is_master())
      {
         iread(5,*n_projector,*nprj);
         iread(5,*l_projector,*nprj);
         iread(5,*m_projector,*nprj);
      }
      parall->Brdcst_iValues(0,0,*nprj,*n_projector);
      parall->Brdcst_iValues(0,0,*nprj,*l_projector);
      parall->Brdcst_iValues(0,0,*nprj,*m_projector);

      nn = (*nmax)*(*nmax)*(*lmax+1);
      *Gijl = new double[nn];
      if (parall->is_master())
      {
         dread(5,*Gijl,nn);
      }
      parall->Brdcst_Values(0,0,nn,*Gijl);
   }
   if (parall->is_master()) dread(5,rcore,1);
   parall->Brdcst_Values(0,0,1,rcore);
   if (*rcore > 0.0)
      *semicore = 1;
   else
      *semicore = 0;


   /* readin vl 3d block */
   tmp2 = new double [mygrid.nfft3d];
   mygrid.t_read(5,tmp2,-1);
   mygrid.t_pack(0,tmp2);
   mygrid.tt_pack_copy(0,tmp2,vl);

   /* reading vnl 3d block */
   if (*nprj > 0) 
   {
      *vnl = new double[(*nprj)*(mygrid.npack(1))];
      prj = *vnl;
      for (i=0; i<(*nprj); ++i)
      {
         mygrid.t_read(5,tmp2,-1);
         mygrid.t_pack(1,tmp2);
         mygrid.tt_pack_copy(1,tmp2,&prj[i*mygrid.npack(1)]);
      }
   }
   if (*semicore)
   {
      nn     = 5*mygrid.npack(0);
      *ncore = new double[nn];
      prj    = *ncore;

      mygrid.t_read(5,tmp2,-1);
      mygrid.t_pack(0,tmp2);
      mygrid.tt_pack_copy(0,tmp2,prj);

      mygrid.t_read(5,tmp2,-1);
      mygrid.t_pack(0,tmp2);
      mygrid.tt_pack_copy(0,tmp2,&prj[2*mygrid.npack(0)]);

      mygrid.t_read(5,tmp2,-1);
      mygrid.t_pack(0,tmp2);
      mygrid.tt_pack_copy(0,tmp2,&prj[3*mygrid.npack(0)]);

      mygrid.t_read(5,tmp2,-1);
      mygrid.t_pack(0,tmp2);
      mygrid.tt_pack_copy(0,tmp2,&prj[4*mygrid.npack(0)]);
   }

   delete [] tmp2;

   if (parall->is_master()) closefile(5);
}
  

/* Constructors */

/*******************************************
 *                                         *
 *     Pseudopotential::Pseudopotential    *
 *                                         *
 *******************************************/
Pseudopotential::Pseudopotential(Ion& myion, PGrid& mygrid)
{
   int ia,version,nfft[3];
   int *n_ptr,*l_ptr,*m_ptr;
   double *rc_ptr,*G_ptr,*vnl_ptr,*ncore_ptr;
   double unita[9];
   char fname[80],aname[2];

   npsp = myion.nkatm;

   psp_type = new int[npsp];
   lmax     = new int[npsp];
   lmmax    = new int[npsp];
   locp     = new int[npsp];
   nmax     = new int[npsp];
   nprj     = new int[npsp];
   semicore = new int[npsp];

   n_projector = new int* [npsp];
   l_projector = new int* [npsp];
   m_projector = new int* [npsp];

   zv          = new double[npsp];
   amass       = new double[npsp];
   rcore       = new double[npsp];
   rc          = new double* [npsp];
   vl          = new double* [npsp];
   for (ia=0; ia<npsp; ++ia) 
      vl[ia] = new double [mygrid.npack(0)];
   Gijl        = new double* [npsp];
   vnl         = new double* [npsp];
   ncore_atom  = new double* [npsp];
   
   comment  = new  char* [npsp];
   for (ia=0; ia<npsp; ++ia) comment[ia] = new char[80];

   for (ia=0; ia<npsp; ++ia)
   {
      strcpy(fname,myion.atom(ia));
      strcat(fname,".vpp");
      psp_read(mygrid,
               fname,
               comment[ia],&psp_type[ia],&version,nfft,unita,aname,
               &amass[ia],&zv[ia],&lmmax[ia],&lmax[ia],&locp[ia],&nmax[ia],
               &rc_ptr,&nprj[ia],&n_ptr,&l_ptr,&m_ptr,&G_ptr,&semicore[ia],&rcore[ia],
               &ncore_ptr,vl[ia],&vnl_ptr);

      rc[ia]          = rc_ptr;
      n_projector[ia] = n_ptr;
      l_projector[ia] = l_ptr;
      m_projector[ia] = m_ptr;
      Gijl[ia]        = G_ptr;
      vnl[ia]         = vnl_ptr;
      if (semicore[ia])
         ncore_atom[ia]  = ncore_ptr;
   }

}

