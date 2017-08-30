

#include	<iostream>
using namespace std;

#include	"Parallel.h"

main(int argc, char *argv[])
{
   double sum,sum0[50];
   int dim,proc_geom[30];
   Parallel parall(argc,argv);

   cout << "np    =" << parall.np_i[0]     << "\n";
   cout << "taskid=" << parall.taskid_i[0] << "\n";
   //cout << "comm  =" << parall.comm   << "\n";
   for (int i=0; i<=parall.dim; ++i)
      cout << "i="<< i << " np_i  =" << parall.np_i[i] << " taskid_i=" << parall.taskid_i[i] << "\n";

   sum = ((double) parall.taskid_i[0]);
   for (int i=0; i<=parall.dim; ++i)
      sum0[i] = parall.SumAll(i,sum);

   cout << "taskid=" << parall.taskid_i[0] << " sum=" << sum;
   for (int i=0; i<=parall.dim; ++i)
        cout << " i=" << i << ",sum0=" << sum0 [i];
   cout << "\n";
}
