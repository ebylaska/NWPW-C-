
#include	<stdint.h>
#include	<stdio.h>
#include	"rtdb_seq.h"

main()
{
   int status;
   int handle;
   int64_t array[30];
   

   printf("sizeof=%d %d\n",sizeof(int),sizeof(long));
   status = rtdb_seq_open("eric.db","old",&handle);
   printf("status=%d\n",status);

   status = rtdb_seq_get(handle,"cpsd:loop",rtdb_int,2,array);
   printf("status2=%d\n",status);
   printf("array=%d %d\n",array[0],array[1]);
}
