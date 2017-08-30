
#include	<iostream>
#include	<cstdio>
#include	<stdio.h>
#include	<cmath>
#include	<cstdlib>
using namespace std;

#include	"Parallel.h"
#include	"control.h"
#include	"util_date.h"
#include	"PGrid.h"
#include	"Ion.h"
#include	"Pseudopotential.h"


int main(int argc, char *argv[])
{
   Parallel myparallel(argc,argv);
   RTDB myrtdb(&myparallel, "eric.db", "old");

   int ii,ia,ngrid[3],matype,nelem;
   char date[26];
   double cpu1,cpu2,cpu3,cpu4;

   if (myparallel.is_master())
   {
      seconds(&cpu1);
      ios_base::sync_with_stdio();
      cout << "          *****************************************************\n";
      cout << "          *                                                   *\n";
      cout << "          *     Car-Parrinello microcluster calculation       *\n";
      cout << "          *                                                   *\n";
      cout << "          *     [     steepest descent minimization   ]       *\n";
      cout << "          *     [          C++ implementation         ]       *\n";
      cout << "          *                                                   *\n";
      cout << "          *            version #6.00   07/20/09               *\n";
      cout << "          *                                                   *\n";
      cout << "          *    This code was developed by Eric J. Bylaska,    *\n";
      cout << "          *    and was based upon algorithms and code         *\n";
      cout << "          *    developed by the group of Prof. John H. Weare  *\n";
      cout << "          *                                                   *\n";
      cout << "          *****************************************************\n";
      cout << "          >>> job started at       " << util_date() << " <<<\n";
   }

   control_read(myrtdb);
   lattice_init();
   myparallel.init2d(control_np_orbital());

   /* initialize lattice, parallel grid structure */
   PGrid mygrid(&myparallel);

  /* initialize psi1 and psi2 */
  //Pneb psi1(&mygrid);
  //Pneb psi2(psi1);

   /* read in ion structure */
   Ion myion(myrtdb);

   Pseudopotential mypsp(myion,mygrid);

   /* setup ewald */
   //Ewald myewald(&myparallel,&myion)



//                 |**************************|
// *****************   summary of input data  **********************
//                 |**************************|

   if (myparallel.is_master())
   {
      cout << "\n";
      cout << " number of processors used: " << myparallel.np() << "\n";
      cout << " processor grid           : " << myparallel.np_i() << " x" << myparallel.np_j() << "\n";
      if (mygrid.maptype==1) cout << " parallel mapping         : slab"    << "\n";
      if (mygrid.maptype==2) cout << " parallel mapping         : hilbert" << "\n";
      if (mygrid.isbalanced()) 
         cout << " parallel mapping         : balanced" << "\n";
      else
         cout << " parallel mapping         : not balanced" << "\n";

      cout << "\n options:\n";
      cout << "\n input movecs:" << control_input_movecs_filename() << "\n";
  
      cout << "\n elements involved in the cluster:\n";
      for (ia=0; ia<myion.nkatm; ++ia)
      {
         printf("    %2d : %4s   core charge: %4.1lf  lmax=%1d\n",
                 ia,myion.atom(ia),mypsp.zv[ia],mypsp.lmax[ia]);
         printf("           comment : %s\n",mypsp.comment[ia]);
         printf("           pseudopotential type            : %3d\n",mypsp.psp_type[ia]);
         printf("           highest angular component       : %3d\n",mypsp.lmax[ia]);
         printf("           local potential used            : %3d\n",mypsp.locp[ia]);
         printf("           number of non-local projections : %3d\n",mypsp.nprj[ia]);
         if (mypsp.semicore[ia])
            printf("           semicore corrections included   : %6.3lf (radius) %6.3lf (charge)\n",mypsp.rcore[ia],mypsp.ncore(ia));
         printf("           cutoff = ");
         for (ii=0; ii<=mypsp.lmax[ia]; ++ii)
            printf("%8.3lf",mypsp.rc[ia][ii]);
         printf("\n");
      }

      cout << "\n atom composition:" << "\n";
      for (ia=0; ia<myion.nkatm; ++ia)
         cout << "   " << myion.atom(ia) << " : " << myion.natm[ia];
      cout << "\n\n initial ion positions (au):" << "\n";
      for (ii=0; ii<myion.nion; ++ii)
         printf("%4d %s\t %8.3lf %8.3lf %8.3lf\n",ii+1,myion.symbol(ii),
                                               myion.rion1[3*ii],
                                               myion.rion1[3*ii+1],
                                               myion.rion1[3*ii+2]);
      cout << "\n";
      cout << " supercell:\n";
      printf("      volume : %10.2lf\n",lattice_omega());
      printf("      lattice:    a1=< %8.3lf %8.3lf %8.3lf >\n",lattice_unita(0,0),lattice_unita(1,0),lattice_unita(2,0));
      printf("                  a2=< %8.3lf %8.3lf %8.3lf >\n",lattice_unita(0,1),lattice_unita(1,1),lattice_unita(2,1));
      printf("                  a3=< %8.3lf %8.3lf %8.3lf >\n",lattice_unita(0,2),lattice_unita(1,2),lattice_unita(2,2));
      printf("      reciprocal: b1=< %8.3lf %8.3lf %8.3lf >\n",lattice_unitg(0,0),lattice_unitg(1,0),lattice_unitg(2,0));
      printf("                  b2=< %8.3lf %8.3lf %8.3lf >\n",lattice_unitg(0,1),lattice_unitg(1,1),lattice_unitg(2,1));
      printf("                  b3=< %8.3lf %8.3lf %8.3lf >\n",lattice_unitg(0,2),lattice_unitg(1,2),lattice_unitg(2,2));

      printf("      density cutoff= %7.3lf fft= %4d x %4d x %4d  (%8d waves %8d per task)\n",
             lattice_ecut(),mygrid.nx,mygrid.ny,mygrid.nz,mygrid.npack_all(0),mygrid.npack(0));
      printf("      wavefnc cutoff= %7.3lf fft= %4d x %4d x %4d  (%8d waves %8d per task)\n",
             lattice_wcut(),mygrid.nx,mygrid.ny,mygrid.nz,mygrid.npack_all(1),mygrid.npack(1));

   }

//                 |**************************|
// *****************   start iterations       **********************
//                 |**************************|
   if (myparallel.is_master()) 
   {
      seconds(&cpu2);
      cout << "          ================ iteration =========================\n";
      cout << "          >>> iteration started at " << util_date() << " <<<\n";

   }
   if (myparallel.is_master()) 
   {
      seconds(&cpu3);
      cout << "          >>> iteration ended at   " << util_date() << " <<<\n";
   }

//                  |***************************|
// ****************** report summary of results **********************
//                  |***************************|
   if (myparallel.is_master()) 
   {
      cout << "\n";
      cout << "          =============  summary of results  =================\n";
      cout << "\n final ion positions (au):" << "\n";
      for (ii=0; ii<myion.nion; ++ii)
         printf("%4d %s\t %8.3lf %8.3lf %8.3lf\n",ii+1,myion.symbol(ii),
                                               myion.rion1[3*ii],
                                               myion.rion1[3*ii+1],
                                               myion.rion1[3*ii+2]);
      cout << "\n\n";
      printf(" total     energy    : %19.10le (%15.5le /ion)\n",      0.0,0.0);
      printf(" total orbtial energy: %19.10le (%15.5le /electron)\n", 0.0,0.0);
      printf(" hartree energy      : %19.10le (%15.5le /electron)\n", 0.0,0.0);
      printf(" exc-corr energy     : %19.10le (%15.5le /electron)\n", 0.0,0.0);
      printf(" ion-ion energy      : %19.10le (%15.5le /ion)\n",      0.0,0.0);

   }

//                 |**************************|
// *****************   report consumed time   **********************
//                 |**************************|
   if (myparallel.is_master()) 
   {
      seconds(&cpu4);
      double t1 = cpu2-cpu1;
      double t2 = cpu3-cpu2;
      double t3 = cpu4-cpu3;
      double t4 = cpu4-cpu1;
      double av = t2/((double ) 100);
      cout << "\n";
      cout << "-----------------"    << "\n";
      cout << "cputime in seconds"   << "\n";
      cout << "prologue    : " << t1 << "\n";
      cout << "main loop   : " << t2 << "\n";
      cout << "epilogue    : " << t3 << "\n";
      cout << "total       : " << t4 << "\n";
      cout << "cputime/step: " << av << "\n";
      cout << "\n";
      cout << "          >>> job completed at     " << util_date() << " <<<\n";

   }
}
