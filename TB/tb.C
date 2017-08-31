/* TB.C -
   Author - Eric Bylaska

*/


#include	<math.h>
#include	"Interactions.h"
#include	"TB.h"


/********************************
 *				*
 *	Constructor		*
 *				*
 ********************************/

inline TB::TB(const int spin, const int nelc, const double U, 
       const int nions_in, Ion* ions_in, Ion* fions_in) 
	: TB_Eigen(4*nions_in, spin, nelc, U)
{
  
   nions = nions_in;
   ions  = ions_in;
   fions = fions_in;
   neigs = 4*nions;

   this->Setup();

}
/********************************
 *				*
 *	     Reset		*
 *				*
 ********************************/
inline void	TB::Reset()
{
   this->Setup();
}

/********************************
 *				*
 *	     Setup		*
 *				*
 ********************************/

inline void TB::Setup()
{
   int 	i,j;
   int  ii,jj; 
   int alpha,beta;
   double tmp;
   

   /* Calculate the Core Energy */
   coreE = 0.0;
   for (i=0; i<nions; ++i)
      for (j=i+1; j<nions; ++j)
         coreE += Ecore(ions[i],ions[j]);

   /* setup the TB Matrix */
   for (i=0;     i<nions;  ++i)
   for (alpha=1; alpha<=4; ++alpha)
   {
      ii = 4*i + (alpha-1);

      /* Onsight elements */
      for (beta=1; beta<=4; ++beta)
      { 
         jj = 4*i + (beta-1);
	 this->Set_Matrix(ii,jj,0.0);
      }
      tmp = Onsight(alpha,ions[i]);
      this->Set_Matrix(ii,ii,tmp);

   
      /* Diagonal Elements */
      for (j=i+1;    j<nions; ++j)
      for (beta=1; beta<=4; ++beta)
      {
         jj = 4*j + (beta-1);
         tmp = Hopping(alpha,ions[i],beta,ions[j]);
         this->Set_Matrix(ii,jj,tmp);
      } /* for j */

   } /* for i*/

   /* solve the TB matrix */
   this->Diagonalize();

   /* set the total energy */
   totalE = coreE 
	  + (this->Valence_Energy()) 
	  + (this->U_Energy());
   
}


/********************************
 *				*
 *	     Force		*
 *				*
 ********************************/

inline void TB::Force()
{
   int i,j;


   for (i=0; i<nions; ++i)
      fions[i] = Vector3(0,0,0);

   /* Core Force */
   for (i=0; i<nions; ++i)
      for (j=i+1; j<nions; ++j)
      {
          Vector3 v     = ions[i] - ions[j];
          double  norm  = 1.0/v.r();

          fions[i] = fions[i] - (norm*dEcoredr(ions[i],ions[j]))*v;
          fions[j] = fions[j] + (norm*dEcoredr(ions[i],ions[j]))*v;
      }


   /* Valence Force */

     for (i=0; i<nions; ++i)
        for (j=i+1; j<nions; ++j)
        {
            Vector3 v     = ions[i] - ions[j];
            double  norm  = 1.0/v.r();
          
            for (int alpha=1; alpha<=4; ++alpha)
            {
               int ni = 4*i + (alpha-1);
  
               for (int beta=1; beta<=4; ++beta)
               {
                  int    nj   = 4*j + (beta-1);
                  double fijr = dHoppingdr(alpha,ions[i],beta,ions[j]);
  
                  for (int nu=0; nu<neigs; ++nu)
                  {
                     double coeff = 2.0
                                  * ((*this)(nj,nu))
                                  * ((*this)(ni,nu))
                                  * (this->Fill(nu));
  	
  
                     fions[i] = fions[i] - v*(norm*coeff*fijr);
                     fions[j] = fions[j] + v*(norm*coeff*fijr);
                  }/*for nu*/ 
               }/*for beta*/
            }/*for alpha*/
         }/*for ij*/


}
