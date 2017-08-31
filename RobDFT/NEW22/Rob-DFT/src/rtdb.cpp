/* rtdb.C
   Author - Eric Bylaska

*/


#include        <iostream>
#include        <cstdlib>
using namespace std;


#include	"Int64.h"
#include	"rtdb.h"

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
RTDB::RTDB(const char *filename, const char *mode)
{
   if (!rtdb_seq_open(filename,mode, &handle))
   {
      cout << "error opening " << filename << " mode=" << mode << "\n";
      return;
   }

}

static int matypesize(const int intype)
{
  switch (intype) {
  case rtdb_char:       /* char */
    return sizeof(char); break;
  case rtdb_int:        /* int */
    return sizeof(Int64); break;
  case rtdb_log:       /* log */
    return sizeof(int); break;
  case rtdb_long:       /* log */
    return sizeof(long); break;
  case rtdb_float:      /* float */
    return sizeof(float); break;
  case rtdb_double:     /* double */
    return sizeof(double); break;
  default:
    return sizeof(char); break;
  }
}


/********************************
 *                              *
 *       RTDB::get              *
 *                              *
 ********************************/
int RTDB::get(const char *tag, const int matype, const int nelem, void *array)
{
   int status;

   status = rtdb_seq_get(handle,tag,matype,nelem,array);
   return status;
}

/********************************
 *                              *
 *       RTDB::put              *
 *                              *
 ********************************/
int RTDB::put(const char *tag, const int matype, const int nelem, void *array)
{
   int status;

   status = rtdb_seq_put(handle,tag,matype,nelem,array);
   if (!status)
      cout << "rtdb error putting tag =" << tag << "\n";
   return status;
}

int RTDB::get_info(const char *tag, int *matype, int *nelem, char *date)
{
   int status;

   status = rtdb_seq_get_info(handle,tag,matype,nelem,date);
   if (!status)
      cout << "rtdb error get_info tag =" << tag << "\n";
   return status;
}
