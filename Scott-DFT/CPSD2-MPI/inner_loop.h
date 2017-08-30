#ifndef _INNERLOOP_H_
#define _INNERLOOP_H_


#include        "Pneb.h"
#include        "PGrid.h"
#include        "Ion.h"
#include        "Ewald.h"
#include	"Kinetic.h"
#include	"Coulomb.h"
#include	"Strfac.h"
#include        "Pseudopotential.h"

extern void inner_loop(Pneb *, Ion *, 
                       Kinetic_Operator *, Coulomb_Operator *,
                       Pseudopotential *, Strfac *, Ewald *,
                       double *, double *, double *, double *,
                       double *, double *, double *,
                       double *, double *, double *, double *);
#endif
