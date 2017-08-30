#ifndef ATLAS_MV_H
   #define ATLAS_MV_H

#include "atlas_misc.h"
#if defined(SREAL)
   #include "atlas_smv.h"
#elif defined(DREAL)
   #include "atlas_dmv.h"
#elif defined(SCPLX)
   #include "atlas_cmv.h"
#elif defined(DCPLX)
   #include "atlas_zmv.h"
#endif

#endif
