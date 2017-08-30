
#include	<iostream>
#include	<cmath>
#include	<cstdlib>
using namespace std;

#include	"olist.h"

main()
{
   OList olist(10);

   olist.insert(5);
   olist.insert(3);
   olist.insert(2);
   olist.print();

   cout << "index 3 = " << olist.index(3) << "\n"; 
   cout << "index 5 = " << olist.index(5) << "\n";
}
