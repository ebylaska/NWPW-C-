/* xyz_io.C
   Author - Eric Bylaska
*/
#include	<iostream.h>
#include	"Ion.h"
#include	"xyz_io.h"

#define	tobohr	0.529177


void	Read_xyz(istream& stream, int* nions, Ion* ions)
{
   int N,i;

   stream >> N;
   *nions = N;


    for (i=0; i<N; ++i)
    {
       stream >> ions[i];
       ions[i] = (1.0/tobohr)*ions[i];
    }
}

void	Write_xyz(ostream& stream, int nions, Ion* ions)
{
   stream << nions << "\n\n";

   for (int i=0; i<nions; ++i)
   {
      ions[i] = tobohr*ions[i];
      stream << ions[i];
      ions[i] = (1.0/tobohr)*ions[i];
   }
}


