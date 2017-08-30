
#include <iostream>
using namespace std;

#include "mpi.h"

int
main(int argc, char *argv[])
{
  MPI::Init(argc, argv);


  int npi,npj,tmp[9],taskid_j;
  
  int rank = MPI::COMM_WORLD.Get_rank();
  int size = MPI::COMM_WORLD.Get_size();
  
  cout << "Hello World! I am " << rank << " of " << size << endl;

  MPI::Intracomm comm2;
  MPI::Group     group2,group3;


  npi = 2; npj = 2;
  if (rank<2) 
    taskid_j=0;
  else 
     taskid_j=1;
  for (int i=0; i<npi; ++i) tmp[i]  = i + 2*taskid_j;

  //group2 = MPI::COMM_WORLD.Get_group();
  //group3 = group2.Incl(2,tmp);
  group2 = MPI::COMM_WORLD.Get_group().Incl(2,tmp);
  comm2  = MPI::COMM_WORLD.Create(group2);

  //worldgroup = worldgroup.Incl(2,tmp);
  int isum2;
  int isum = rank;

  comm2.Allreduce(&isum,&isum2,1,MPI_INTEGER,MPI_SUM);
  
  cout << "isum=" << isum << " isum2=" << isum2 << " taskid_i=" << comm2.Get_rank() << " npi=" << comm2.Get_size() << "\n";
  MPI::Finalize();
}


 

