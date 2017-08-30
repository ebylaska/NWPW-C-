
#include	<iostream>
#include	<cmath>
#include	<cstdlib>
using namespace std;


#include	"Parallel.h"
#include	"d3db.h"
#include	"compressed_io.h"

main(int argc,char *argv[])
{
   int j[100];
   double *psi;
   int taskid,np,nx,ny,nz;
   np = 3; taskid=0;
   nx = ny = nz = 32;

   Parallel parall(argc,argv);
   d3db map4(parall,2,nx,ny,nz);
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

   openfile(6,"DATAFILE","w",strlen("DATAFILE"));
      for (int i=0; i<10; ++i) j[i] = i*10;
      iwrite(6,j,10);
   closefile(6);

   openfile(5,"DATAFILE","r",strlen("DATAFILE"));
   for (int i=0; i<10; ++i)
   {
      iread(5,j,1);
      cout << i << ",j=" << j[0] << "\n";
   }
   closefile(5);




}
