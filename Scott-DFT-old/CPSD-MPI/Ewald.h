#ifndef	_EWALD_H_
#define _EWALD_H_

using namespace std;

#include        "Parallel.h"
#include        "lattice.h"
#include	"Ion.h"

class	Ewald {
   int    ncut,enx,eny,nyz,nshl3d,enpack,nida;
   int    *i_indx,*j_indx,*k_indx;
   double *vg,*rcell,*eG,*vcx,*zv;
   double *ewx1,*ewx2,*ewx3;
   double rcut,cewald,alpha;

public:
   Parallel  *ewaldparall;
   Ion	     *ewaldion;

   /* Constructors */
   Ewald(Parallel *, Ion *);

   /* destructor */
   ~Ewald() {
            delete [] i_indx;
            delete [] j_indx;
            delete [] k_indx;
            delete [] vg;
            delete [] rcell;
            delete [] eG;
            delete [] zv;
            delete [] ewx1;
            delete [] ewx2;
            delete [] ewx3;
         }

    }


    void phafac();
    int ncut();
    int nida()
    int npack();
    int nshl3d()
    double zv();
    double rcut();
    double mandelung();
    double e();
    void   f(double *);

};

#endif
