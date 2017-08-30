#ifndef ATLAS_NTHREADS_H
   #define ATLAS_NTHREADS_H

/* Get rid of 00 if you don't want to build pthreads */
   #ifndef ATL_WINTHREADS00
      #include "pthread.h"
   #endif
   #define ATL_NTHREADS 8
   #define ATL_NTHRPOW2 3
   #ifdef ATL_LAUNCHORDER
       static int ATL_launchorder[8] = {0,4,2,6,1,5,3,7};
   #endif

#endif
