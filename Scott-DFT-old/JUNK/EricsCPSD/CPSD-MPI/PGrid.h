#ifndef _PGRID_H_
#define _PGRID_H_
/* d3db.h
   Author - Eric Bylaska

*/

#include	"Parallel.h"
#include	"d3db.h"
#include	"lattice.h"
#include	"Balance.h"

class PGrid : public d3db, public lattice {

   Balance *mybalance;
   int balanced;
   double *Garray;
   int     *masker[2],*packarray[2];
   int     nwave[2],nwave_entire[2],nwave_all[2],nida[2],nidb[2],nidb2[2];

public:

        /* constructor */
	PGrid(Parallel *);

        /* destructor */
        ~PGrid() { 
            delete [] Garray; 
            delete [] masker[0];
            delete [] packarray[0];
            delete mybalance;
         }

         int npack(const int i)     {return (nida[i] + nidb[i]);}
         int npack_all(const int i) {return nwave_all[i];}
         int isbalanced() { return balanced;  }
         void c_pack(const int, double *a);
};

#endif
