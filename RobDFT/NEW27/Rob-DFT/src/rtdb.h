#ifndef _RTDB_H_
#define _RTDB_H_
/* rtdb.h
   Author - Eric Bylaska

	this class is used defining nd parallel geometries
*/

#include	"rtdb_seq.h"

class RTDB {

   int handle;

public:

	/* Constructors */
	RTDB(const char *, const char *);

        /* destructor */
	//~RTDB();

        int get(const char *, const int, const int, void *);
        int put(const char *, const int, const int, void *);
        int get_info(const char *, int *, int *, char *);

};
#endif
