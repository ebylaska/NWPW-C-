#ifndef _PNEB_H_
#define _PNEB_H_
/* Pneb.h
   Author - Eric Bylaska

*/

#include	"Parallel.h"
#include	"d3db.h"
#include	"lattice.h"
#include	"Balance.h"
#include	"Pgrid.h"

class Pneb : public PGrid, public d3db, public lattice {

   int ispin,ne[2],neq[2];

public:

        /* constructor */
	Pneb(PGrid *);

        /* destructor */
        ~Pneb() { 
         }

};

#endif
