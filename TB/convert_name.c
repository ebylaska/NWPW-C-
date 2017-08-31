/* convert_name.c - 
   Author - Eric Bylaska

*/
#include	<string.h>

/********************************
 *				*
 *	NameToCharge		*
 *				*
 ********************************/

int	NameToCharge(char* name)
{
   int charge = 0;

   if (strcmp("H",name)  == 0)  charge = 1;
   if (strcmp("He",name) == 0)  charge = 2;
   if (strcmp("Li",name) == 0)  charge = 3;
   if (strcmp("Be",name) == 0)  charge = 4;
   if (strcmp("B",name)  == 0)  charge = 5;
   if (strcmp("C",name)  == 0)  charge = 6;
   if (strcmp("N",name)  == 0)  charge = 7;
   if (strcmp("O",name)  == 0)  charge = 8;
   if (strcmp("F",name)  == 0)  charge = 9;
   if (strcmp("Ne",name) == 0)  charge = 10;
   if (strcmp("Na",name) == 0)  charge = 11;
   if (strcmp("Mg",name) == 0)  charge = 12;
   if (strcmp("Al",name) == 0)  charge = 13;
   if (strcmp("Si",name) == 0)  charge = 14;
   if (strcmp("P",name)  == 0)  charge = 15;
   if (strcmp("S",name)  == 0)  charge = 16;
   if (strcmp("Cl",name) == 0)  charge = 17;
   if (strcmp("Ar",name) == 0)  charge = 18;

   return charge;
}

/********************************
 *				*
 *	NameToCoreCharge	*
 *				*
 ********************************/

int	NameToCoreCharge(char* name)
{
   int charge = 0;

   if (strcmp("H",name)  == 0)  charge = 0;
   if (strcmp("He",name) == 0)  charge = 0;
   if (strcmp("Li",name) == 0)  charge = 2;
   if (strcmp("Be",name) == 0)  charge = 2;
   if (strcmp("B",name)  == 0)  charge = 2;
   if (strcmp("C",name)  == 0)  charge = 2;
   if (strcmp("N",name)  == 0)  charge = 2;
   if (strcmp("O",name)  == 0)  charge = 2;
   if (strcmp("F",name)  == 0)  charge = 2;
   if (strcmp("Ne",name) == 0)  charge = 2;
   if (strcmp("Na",name) == 0)  charge = 6;
   if (strcmp("Mg",name) == 0)  charge = 6;
   if (strcmp("Al",name) == 0)  charge = 6;
   if (strcmp("Si",name) == 0)  charge = 6;
   if (strcmp("P",name)  == 0)  charge = 6;
   if (strcmp("S",name)  == 0)  charge = 6;
   if (strcmp("Cl",name) == 0)  charge = 6;
   if (strcmp("Ar",name) == 0)  charge = 6;

   
   return charge;
}


/********************************
 *				*
 *	NameToMass  		*
 *				*
 ********************************/

double	NameToMass(char* name)
{
   double mass = 0.0;

   if (strcmp("H",name)  == 0)  mass   = 1.0;
   if (strcmp("He",name) == 0)  mass   = 4.0;
   if (strcmp("Li",name) == 0)  mass   = 7.0;
   if (strcmp("Be",name) == 0)  mass   = 9.0;
   if (strcmp("B",name)  == 0)  mass   = 11.0;
   if (strcmp("C",name)  == 0)  mass   = 12.0;
   if (strcmp("N",name)  == 0)  mass   = 14.0;
   if (strcmp("O",name)  == 0)  mass   = 16.0;
   if (strcmp("F",name)  == 0)  mass   = 19.0;
   if (strcmp("Ne",name) == 0)  mass   = 20.0;
   if (strcmp("Na",name) == 0)  mass   = 23.0;
   if (strcmp("Mg",name) == 0)  mass   = 24.0;
   if (strcmp("Al",name) == 0)  mass   = 27.0;
   if (strcmp("Si",name) == 0)  mass   = 28.0;
   if (strcmp("P",name)  == 0)  mass   = 31.0;
   if (strcmp("S",name)  == 0)  mass   = 32.0;
   if (strcmp("Cl",name) == 0)  mass   = 35.0;
   if (strcmp("Ar",name) == 0)  mass   = 40.0;


   return mass;
}


/********************************
 *				*
 *	NameToValenceCharge	*
 *				*
 ********************************/

int	NameToValenceCharge(char* name)
{
   int q = NameToCharge(name) - NameToCoreCharge(name);

   return q;
}


