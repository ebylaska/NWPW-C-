
using namespace std;

#include	<fstream.h>
#include	"error.h"
#include	"Slater_Koster.h"
#include	"Slater_Koster_List.h"
#include	"Projection.h"
#include	"Interactions.h"

static Slater_Koster_List	list;

void Read_Interactions(ifstream& stream)
{
    stream >> list;
}
double	Onsight(const int alpha, Ion& ion)
{
   double tmp = Interaction(ion,ion)->t_aa(alpha);

   return tmp;
}

double	Hopping(const int alpha, Ion& ion_a,
		const int beta,  Ion& ion_b)
{
   Projection	p(ion_a,ion_b);

   double	tmp = Interaction(ion_a,ion_b)->t_ab(alpha,beta,p);

   return tmp;
}

double	dHoppingdr(const int alpha, Ion& ion_a,
		const int beta,  Ion& ion_b)
{
   Projection	p(ion_a,ion_b);
   
   double	tmp = Interaction(ion_a,ion_b)->dt_abdr(alpha,beta,p);

   return tmp;
}

double	Ecore(Ion& ion_a, Ion& ion_b)
{
   Vector3 v1  = ion_a - ion_b;
   double  r12 = v1.r();
   double  tmp =  Interaction(ion_a,ion_b)->Ecore(r12);

   return tmp;
}

double	dEcoredr(Ion& ion_a, Ion& ion_b)
{
   Vector3 v1  = ion_a - ion_b;
   double  r12 = v1.r();
   double  tmp =  Interaction(ion_a,ion_b)->dEcoredr(r12);

   return tmp;
}



Slater_Koster	*Interaction(Ion& ion1, Ion& ion2)
{
   int			i,N;
   int			zz1 = ion1.Charge();
   int			zz2 = ion2.Charge();
   int			z1,z2;
   int			error_here, found;
   Slater_Koster	*sk;
   

   /* make sure z1 and z2 are in the right order */
   if (zz1 <= zz2) {  z1 = zz1; z2 = zz2;}
   else            {  z1 = zz2; z2 = zz1;}

   N = list.Size();

   i = 0;
   found      = 0;
   error_here = 0;
   while ((!found) && (!error_here))
   {
      sk = list(i);
      if ( (z1 == (sk->Z1())) && (z2 == (sk->Z2())) )
         found = 1;

      ++i;
      if (i >= N) error_here = 1;
   }

      
   /* Error */
   if ((error_here) && (!found))
   {
      char err_string[80];
      strcpy(err_string,"Interaction: interaction(");
      strcat(err_string,ion1.Name());
      strcat(err_string," ");
      strcat(err_string,ion2.Name());
      strcat(err_string,") not found");
      //Error_Found(err_string);
      printf("z1=%d z2=%d  %d %d\n",z1,z2,sk->Z1(),sk->Z2());
   }


    return sk;
}

int	Number_Interactions()
{
    return list.Size();
}

void	Print_Interactions()
{
   Slater_Koster	*sk;
   sk = list(0);

   printf("Hello world\n");
   printf("list size=%d\n",list.Size());
   printf("Z1=%d Z2=%d\n",sk->Z1(),sk->Z2());
   printf("A1=%s A2=%s\n",sk->Atom1(),sk->Atom2());
   printf("F0=%lf %lf %lf %lf %lf \n",sk->F0(0),sk->F0(1),sk->F0(2),sk->F0(3),sk->F0(4));



}
