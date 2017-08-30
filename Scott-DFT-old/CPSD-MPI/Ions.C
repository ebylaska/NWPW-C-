/* Ions.C - 
   Author - Eric Bylaska
*/

using namespace std;

#include	<string.h>
#include	"Convert_Name.h"
#include	"Ions.h"


/* Constructors */

Ion::Ion() 
{
   strcpy(name,"X");
   charge = NameToCharge(name);
   mass   = NameToMass(name);
}

