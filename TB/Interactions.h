#ifndef _Interactions_H_
#define _Interactions_H_

using namespace std;

#include	<fstream.h>
#include	"Slater_Koster.h"
#include	"Ion.h"


extern void 		Read_Interactions(ifstream&);
extern void		Print_Interactions();
extern Slater_Koster	*Interaction(Ion&, Ion&);
extern double		Onsight(const int, Ion&);
extern double		Hopping(const int, Ion&, 
       				const int, Ion&);
extern double		dHoppingdr(const int, Ion&, 
       				   const int, Ion&);
extern double		Ecore(Ion&,  Ion&);
extern double		dEcoredr(Ion&,  Ion&);
extern int		Number_Interactions();

#endif
