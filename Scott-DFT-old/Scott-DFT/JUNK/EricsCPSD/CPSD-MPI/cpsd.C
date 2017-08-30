
#include	<iostream>
#include	<cstdio>
#include	<stdio.h>
#include	<cmath>
#include	<cstdlib>
using namespace std;

#include	"Parallel.h"
#include	"rtdb.h"
#include	"control.h"
#include	"util_date.h"
#include	"d3db.h"
#include	"lattice.h"
#include	"PGrid.h"



int main(int argc, char *argv[])
{
   Parallel myparallel(argc,argv);
   RTDB myrtdb(&myparallel, "eric.db", "old");
   int ngrid[3],matype,nelem;
   char date[26];
   double cpu1,cpu2,cpu3,cpu4;
   if (myparallel.is_master()) seconds(&cpu1);
   

   if (myparallel.is_master())
   {
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
   myparallel.init2d(control_np_orbital());

   /* initialize lattice, parallel grid structure */
   PGrid mygrid(&myparallel);

   //d3db myd3db(&myparallel,control_mapping(),control_ngrid(0),control_ngrid(1),control_ngrid(2));
   //lattice mylattice;


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

      cout << "\n";
      cout << " supercell:\n";
      printf("      volume : %10.2lf\n",mygrid.omega());
      printf("      lattice:    a1=< %8.3lf %8.3lf %8.3lf >\n",mygrid.unita(0,0),mygrid.unita(1,0),mygrid.unita(2,0));
      printf("                  a2=< %8.3lf %8.3lf %8.3lf >\n",mygrid.unita(0,1),mygrid.unita(1,1),mygrid.unita(2,1));
      printf("                  a3=< %8.3lf %8.3lf %8.3lf >\n",mygrid.unita(0,2),mygrid.unita(1,2),mygrid.unita(2,2));
      printf("      reciprocal: b1=< %8.3lf %8.3lf %8.3lf >\n",mygrid.unitg(0,0),mygrid.unitg(1,0),mygrid.unitg(2,0));
      printf("                  b2=< %8.3lf %8.3lf %8.3lf >\n",mygrid.unitg(0,1),mygrid.unitg(1,1),mygrid.unitg(2,1));
      printf("                  b3=< %8.3lf %8.3lf %8.3lf >\n",mygrid.unitg(0,2),mygrid.unitg(1,2),mygrid.unitg(2,2));

      printf("      density cutoff= %7.3lf fft= %4d x %4d x %4d  (%8d waves %8d per task)\n",
             mygrid.ecut(),mygrid.nx,mygrid.ny,mygrid.nz,mygrid.npack_all(0),mygrid.npack(0));
      printf("      wavefnc cutoff= %7.3lf fft= %4d x %4d x %4d  (%8d waves %8d per task)\n",
             mygrid.wcut(),mygrid.nx,mygrid.ny,mygrid.nz,mygrid.npack_all(1),mygrid.npack(1));

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
      cout << "          =============  summary of results  =================\n";
      cout << " final position of ions:\n";
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

   return(0);
}



