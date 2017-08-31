/* TB_Eigen.C
   Author - Eric Bylaska

	This class is used for setting up and solving
a real symmetric eigenvalue problem.

*/
#include	"error.h"
#include	"TB_Eigen.h"

#include	"eigen.h"



/********************************
 *				*
 *        Constructors          *
 *				*
 ********************************/

inline TB_Eigen::TB_Eigen()
{
   N    = 0;
   ispin = 0;
   nelc = 0;
   U    = 0.0;
}

inline TB_Eigen::TB_Eigen(const int tN, 
	           const int tspin, const int tnelc,
	           const double tU)
{
   ispin  = tspin;
   nelc   = tnelc;
   U      = tU;

   N      = tN;
   fill   = new int[N];
   eigens = new double[N];
   H      = new double[N*N];
   psi    = H;
   
}

/********************************
 *				*
 *	  Destructor		*
 *				*
 ********************************/

inline TB_Eigen::~TB_Eigen()
{
   delete [] fill;
   delete [] H;
   delete [] eigens;
}

	
/* routines for setting up the matrix and viewing it */

/********************************
 *				*
 *	    Set_Matrix		*
 *				*
 ********************************/

inline void	TB_Eigen::Set_Matrix(const int i, const int j, 
				const double value)
{
   H[i+j*N] = value; /* column ordered, like fortran */
   H[j+i*N] = value; /* symmetric matrix             */
}

/********************************
 *				*
 *           () operators 	*
 *				*
 ********************************/

inline double	TB_Eigen::operator () (const int i, const int j)
{
   return( H[i+j*N]);
}

inline double	TB_Eigen::operator () (const int k)
{
   return( eigens[k] );
}

/********************************
 *				*
 *        Equality operator 	*
 *				*
 ********************************/

inline TB_Eigen& TB_Eigen::operator =(TB_Eigen& source)
{
   int i;

   if (this != &source)
   {
      delete [] fill;
      delete [] H;
      delete [] eigens;

      ispin  = source.ispin;
      spin   = source.spin;
      nelc   = source.nelc;
      U      = source.U;
      Evalence = source.Evalence;
      Uterm    = source.Uterm;

      N      = source.N;
      fill   = new int[N];
      eigens = new double[N];
      H      = new double[N*N];
      psi    = H;
      
      for (i=0; i<N*N; ++i) H[i]      = source.H[i];
      for (i=0; i<N; ++i)   eigens[i] = source.eigens[i];
      for (i=0; i<N; ++i)   fill[i]   = source.fill[i];
   }

   return *this;
}


/********************************
 *				*
 *	   Diagonalize		*
 *				*
 ********************************/

inline void	TB_Eigen::Diagonalize()
{
   int i;

   /* Diagonalize the Matrix */

   int    info  = 0;
   double * work = new double[3*N];
   EIGSOLVER(&N, &N,H,eigens,work,&info);
   delete [] work;

   /* Error in eigenvalue solver */
   if (info != 0)
      Error_Found("EigenSolve::Diagonalize: error in solver");


   /* Figure out filling and energies */
   for (i=0; i<N; ++i)
      fill[i] = 0;

   if (ispin==0)
   {
      int up_counter   = N-1;
      int down_counter = N-1;
      for (i=0; i<nelc; ++i)
      {
         double e1 = eigens[up_counter];
         double e2 = eigens[down_counter] + U;
         if (e1 <= e2)
         {
            fill[up_counter] = fill[up_counter] + 1;
            up_counter       = up_counter - 1;
         }
         else
         {
            fill[down_counter] = fill[down_counter] + 1;
            down_counter       = down_counter - 1;
         }
      }
      spin = (down_counter - up_counter) + 1;
   } /*if*/
   /* ispin!= 0 */
   else
   {
      spin = ispin;
      int nup   = (spin-1+nelc)/2;
      int ndown = nelc-nup;
      int bot1  = N-1;
      int bot2  = N-1;
      for (i=0; i<nup; ++i)
      {
         fill[bot1] = fill[bot1] + 1;
         bot1       = bot1 - 1;
      }
      for (i=0; i<ndown; ++i)
      {
         fill[bot2] = fill[bot2] + 1;
         bot2       = bot2 - 1;
      }
   } /*else*/
     

   /* add up the eigenvalues */
   Evalence = 0.0;
   Uterm    = 0.0;
   for (i=0; i<N; ++i)
   {
      Evalence += fill[i]*eigens[i];
      if (fill[i] > 1)
      Uterm    += U;
   }
   

} /* Diagonalize */

inline int	TB_Eigen::Fill(const int i)
{
   return fill[i];
}

inline double	TB_Eigen::Valence_Energy()
{
   return Evalence;
}
inline double	TB_Eigen::U_Energy()
{
   return Uterm;
}







/********************************
 *				*
 *       operator <<		*
 *				* 
 ********************************/

inline void	operator <<(ostream& s, TB_Eigen& source)
{
   int i,j;

   /* first put out eigenvalues */
   s << "**********************************************\n";
   s << "EigenValues: \n";
   for (i=0; i<source.N; ++i)
      s << source.eigens[i] << "\n";
   s << "\n";

   /* Then Matrix */
   s << source.N << "x" << source.N  << " EigenMatrix: \n";
   for (i=0; i<source.N; ++i)
   {
      for (j=0; j<source.N; ++j)
         s << source.H[i+j*(source.N)] << " ";
      s << "\n";
   }
   s << "**********************************************\n";
}
       


