#ifndef	_XYZ_IO_
#define _XYZ_IO_
/* xyz_io.h
   Author - Eric Bylaska
*/
#include	<iostream.h>
#include	"Ion.h"


void	Read_xyz(istream&,  int*, Ion*);
void	Write_xyz(ostream&, int, Ion*);

#endif

