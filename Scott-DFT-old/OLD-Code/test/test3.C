
#include	<iostream>
#include	<cmath>
#include	<cstdlib>
using namespace std;

#include	"Parallel.h"

main(int argc, char *argv[])
{
   int dim,proc_geom[30];
   Parallel parall(argc,argv);

   cout << "np    =" << parall.np     << "\n";
   cout << "taskid=" << parall.taskid << "\n";
   cout << "comm  =" << parall.comm   << "\n";
   for (int i=0; i<parall.dim; ++i)
      cout << "i="<< i << " np_i  =" << parall.np_i[i] << "\n";

}
