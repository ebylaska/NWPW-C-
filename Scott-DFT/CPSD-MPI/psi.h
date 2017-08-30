#ifndef _PSIGETHEADER_H_
#define _PSIGETHEADER_H_

#include	"compressed_io.h"
#include	"control.h"

#include	"Parallel.h"
#include	"Pneb.h"

extern void psi_get_header(Parallel *, int *, int *, double *, int *, int *);
extern void psi_read(Pneb *, int *, int *, double *, int *, int *,double *);

#endif
