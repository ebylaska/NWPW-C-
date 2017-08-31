#ifndef	_PSEUDOPOTENTIAL_H_
#define _PSEUDOPOTENTIAL_H_


using namespace std;

#include	"Ion.h"
#include	"Pneb.h"
#include	"Strfac.h"

class	Pseudopotential {

   int nprj_max;

   double **Gijl;
   double **ncore_atom;
   double **vl;
   double **vnl;
   double *rlocal;
   double *amass;

   //char **atomsym;
   Pneb   *mypneb;
   Ion    *myion;
   Strfac *mystrfac;

public:
   int npsp;
   int *nprj,*lmax,*lmmax,*locp,*nmax,*psp_type,*semicore;
   int **n_projector,**l_projector,**m_projector;
   double **rc;
   double *zv,*rcore;
   char **comment;

   /* Constructors */
   Pseudopotential(Ion *, Pneb *, Strfac *);

   /* destructor */
   ~Pseudopotential() {}

    double ncore(const int ia) {return 0.0;}

    void v_nonlocal(double *, double *);
    void v_local(double *);

};

#endif
