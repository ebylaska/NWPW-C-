/*
 $Id: compressed_io.c,v 1.3 2004-05-05 02:20:33 edo Exp $
*************************************************************************
*									*
* compressed_io.c							*
*									*
* These routines implement a simple input/output package for writing to	*
* and reading from gzipped files.					*
*                                                                       *
*									*
* Fortran usage:							*
*	call openfile(unit,'file', 'X', n)	X is either 'r' or 'w'	*
*					n is the length of the filename	*
*	call cwrite(unit,c, n)		c is n element character array	*
*	call iwrite(unit,i, n)		i is n element integer array	*
*	call dwrite(unit,d, n)		d is n element double array	*
*	call cread(unit,c, n)		c is n element character array	*
*	call iread(unit,i, n)		i is n element integer array	*
*	call dread(unit,d, n)		d is n element double array	*
*	call closefile(unit)		close the datafile		*
*									*
* Author:      Scott Kohn (skohn@chem.ucsd.edu)				*
* modified by: Eric Bylaska (ebylaska@chem.ucsd.edu)			*
*									*
*************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>



#define MAX_UNIT	10

static FILE* fd[MAX_UNIT];	/* the file descriptor of the pipe */

#define BAIL(X) { fprintf(stderr, X); exit(-1); }

/*
*************************************************************************
*									*
* Define the Xwrite and Xread routines using the Fortran bindings.	*
*									*
*************************************************************************
*/

void cwrite(int unit, char *c, const int n)
{

   (void) fwrite(c, sizeof(char), n, fd[unit]);
}

void cread(int unit, char *c, const int n)
{
   (void) fread(c, sizeof(char), n, fd[unit]);
}

void iwrite(const int unit, const int *i, const int n)
{
   (void) fwrite(i, sizeof(int), n, fd[unit]);
}

void iread(const int unit, int *i, const int n)
{
   (void) fread(i, sizeof(int), n, fd[unit]);
}

void dwrite(const int unit, const double *d, const int n)
{
   (void) fwrite(d, sizeof(double), n, fd[unit]);
}

void dread(const int unit, double *d, const int n)
{
   (void) fread(d, sizeof(double), n, fd[unit]);
}

/*
*************************************************************************
*									*
* void openfile(char *filename, char *mode, Integer *n)			*
* void closefile()							*
*									*
* Function openfile opens a pipe to either gzip (to compress a stream)	*
* or zcat (to uncompress a stream).  Function closefile() closes the	*
* pipe stream created by openfile().					*
*									*
*************************************************************************
*/

#define FUDGE_FACTOR (8)

void openfile(const int unit, char *filename, char *mode, const int n1)
{
   char *file = (char *) malloc(n1+1);
   (void) strncpy(file, filename, n1);
   file[n1] = '\0';

   if ((*mode == 'r') || (*mode == 'R')) {
      if (!(fd[unit] = fopen(file, "rb")))
         BAIL("ERROR:  Could not open pipe from input file\n");
   } else {
      if (!(fd[unit] = fopen(file, "wb")))
         BAIL("ERROR:  Could not open pipe to output file\n");
   }
   free(file);
}

void closefile(const int unit)
{
   (void) fclose(fd[unit]);
}
