/* seconds.C -
   Author - Eric Bylaska

*/
#include	<iostream.h>
#include	<time.h>
#include	"seconds.h"

void	seconds(double* tt)
{
   *tt = Clock/(CLOCKS_PER_SEC);

}
