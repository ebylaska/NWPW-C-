
#include	<iostream>
#include	<cmath>
#include	<cstdlib>
using namespace std;

#include	"Parallel.h"

main(int argc, char *argv[])
{
   int dim,proc_geom[30];
   proc_geom[0] = -1; proc_geom[1] = -1; proc_geom[2] = -1;

   Parallel parall(argc,argv,dim,proc_geom);

   cout << "np    =" << parall.np     << "\n";
   cout << "taskid=" << parall.taskid << "\n";
   cout << "comm  =" << parall.comm   << "\n";
   cout << "np_i  =" << parall.np_i[0] << "\n";
   cout << "np_j  =" << parall.np_i[1] << "\n";

}
