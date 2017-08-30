
#include	<iostream>
#include	<cmath>
#include	<cstdlib>
using namespace std;

#include	<mpi.h>

#include	"d3db.h"
#include	"Mapping3.h"

main(int argc,char *argv[])
{
   double *psi;
   int taskid,np,nx,ny,nz;
   np = 3; taskid=0;
   nx = ny = nz = 32;

   MPI_Init(&argc,&argv);
   d3db map4(2,MPI_COMM_WORLD,nx,ny,nz);
   psi = map4.r_alloc();
   map4.r_zero(psi);


   if (map4.taskid==0) 
   {
   cout << "\n\n";
   cout << "map4type=" << map4.maptype << "\n";
   cout << "nfft3d  =" << map4.nfft3d  << "\n";
   cout << "n2ft3d  =" << map4.n2ft3d  << "\n";
   cout << "nx      =" << map4.nx  << "\n";
   cout << "ny      =" << map4.ny  << "\n";
   cout << "nz      =" << map4.nz  << "\n";
   cout << "nq1     =" << map4.nq1  << "\n";
   cout << "nq2     =" << map4.nq2  << "\n";
   cout << "nq3     =" << map4.nq3  << "\n";
   cout << "np      =" << map4.np   << "\n"; 
   }

   map4.r_dealloc(psi);

   MPI_Finalize();
}
