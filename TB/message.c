/* message.c -
   Author - Eric Bylaska

*/
#include	<stdio.h>
#include	<string.h>
#include	"convert_name.h"
#include	"message.h"

static	char* header1 =
"             ********************************************\n";
static	char* header2 =
"             *                                          *\n";
static	char* header3 =
"             *           Tight Binding Code             *\n";
static	char* header4 =
"             *                                          *\n";
static	char* header5 =
"             *           (steepest descent)             *\n";
static	char* header6 =
"             *             (version 1.0)                *\n";
static	char* header7 =
"             *                                          *\n";
static	char* header8 =
"             ********************************************\n\n";

void	Header_Message()
{
   printf("%s",header1);
   printf("%s",header2);
   printf("%s",header3);
   printf("%s",header4);
   printf("%s",header5);
   printf("%s",header6);
   printf("%s",header7);
   printf("%s",header8);
}

static	char* header11 =
"             ********************************************\n";
static	char* header12 =
"             *                                          *\n";
static	char* header13 =
"             *           Tight Binding Code             *\n";
static	char* header14 =
"             *                                          *\n";
static	char* header15 =
"             *          (molecular dynamics)            *\n";
static	char* header16 =
"             *             (version 1.0)                *\n";
static	char* header17 =
"             *                                          *\n";
static	char* header18 =
"             ********************************************\n\n";

void	Header_TBMD_Message()
{
   printf("%s",header11);
   printf("%s",header12);
   printf("%s",header13);
   printf("%s",header14);
   printf("%s",header15);
   printf("%s",header16);
   printf("%s",header17);
   printf("%s",header18);
}



void	Input_Data_Message(const int nelc,  const int ispin,
			   const double dt, const double tole,
			   const int nions, char** inames, double** ions,
			   const char* parameterization)
{
printf(
 "        |-\\____|\\/-----\\/\\/->    Input Data    <-\\/\\/-----\\/|____/-|\n"
);
printf("\n");
printf("  Tight Binding Parameterization: %s\n", parameterization);
printf("\n");

   ELEMENTS_Message(nions,inames,ions);
   ATOMS_Message(nions,inames,ions);
   ELC_Message(nelc,ispin,nions,inames,ions);
   PARAMS_Message(dt,tole);
}

void	Input_Data_TBMD_Message(const int nelc,  const int ispin,
			   const double dt, const double scalev,
			   const int nions, char** inames, 
			   double** ions, double** vions,
 			   const double ke0,      const double ke0_cm,
			   const double ke_start, const double ke_start_cm,
			   const char* parameterization)
{
printf(
 "        |-\\____|\\/-----\\/\\/->    Input Data    <-\\/\\/-----\\/|____/-|\n"
);
printf("\n");
printf("  Tight Binding Parameterization: %s\n", parameterization);
printf("\n");

   ELEMENTS_Message(nions,inames,ions);
   ATOMS_Message(nions,inames,ions);
   VELCS_Message(nions,inames,vions);
   printf("  initial kinetic energy: ion= %le c.o.m.= %le\n",ke0,ke0_cm); 
   printf("  after scaling velocity: ion= %le c.o.m.= %le\n",ke_start,ke_start_cm); 
   printf("  increased energy      : ion= %le c.o.m.= %le\n",
		ke_start-ke0,ke_start_cm-ke0_cm); 
   printf("\n");
   ELC_Message(nelc,ispin,nions,inames,ions);
   PARAMS_TBMD_Message(dt,scalev);
}




void	Output_Data_Message(const double eall, 
                            const double ec, 
                            const double ev,
                            const double eu,
			    const int nions, char** inames, double** ions)
			   
{
printf(
"        |-\\____|\\/-----\\/\\/->   Output Data    <-\\/\\/-----\\/|____/-|\n"
);
printf("\n");
ATOMS_Message(nions,inames,ions);
printf("\n");
printf("     Total Energy  : %14.6le\n",eall);
printf("     Core Energy   : %14.6le\n",ec);
printf("     Valence Energy: %14.6le\n",ev);
printf("     U term Energy : %14.6le\n",eu);
printf("\n");


}

void	Output_Data_TBMD_Message(const double eall, 
                            const double ec, 
                            const double ev,
                            const double eu,
			    const int nions, char** inames, 
			    double** ions, double** vions)
			   
{
printf(
"        |-\\____|\\/-----\\/\\/->   Output Data    <-\\/\\/-----\\/|____/-|\n"
);
printf("\n");
ATOMS_Message(nions,inames,ions);
VELCS_Message(nions,inames,vions);
printf("\n");
printf("     Total Energy  : %14.6le\n",eall);
printf("     Core Energy   : %14.6le\n",ec);
printf("     Valence Energy: %14.6le\n",ev);
printf("     U term Energy : %14.6le\n",eu);
printf("\n");


}

void	Start_Iteration_Message()
{
printf(
"        |-\\____|\\/-----\\/\\/-> iteration started<-\\/\\/-----\\/|____/-|\n"
);
}

void	End_Iteration_Message(const int option)
{
   if (option)
printf(
"        |-\\____|\\/-----\\/\\/-> tolerences ok    <-\\/\\/-----\\/|____/-|\n"
);
   else

printf(
"        |-\\____|\\/-----\\/\\/-> iteration ended  <-\\/\\/-----\\/|____/-|\n"
);

printf("\n\n");
}

/* this string used by the routines following */
static	char dummy_string[500];
static	char dummy2_string[50];

void	Iteration_Message(const int step, const int n, const double* e)
{
   int i;
   printf("%6d   ",step);
   for (i=0; i<n; ++i)
   {
      printf("%13.6le  ",e[i]);
   }
   printf("\n");
}

void	Start_Eigenvalue_Message(const int spin)
{
   printf( " orbital energies: 2*S+1=%d\n", spin);
   printf( " ------------------------------\n");
}

void	Eigenvalue_Message(const int state,const int fill, const double e)
{
   if (fill == 0)
      printf("%3d ****       ****  %12.6lf\n",state,e);
   else if (fill == 1)
      printf("%3d ****   +   ****  %12.6lf\n",state,e);
   else 
      printf("%3d ****  + -  ****  %12.6lf\n",state,e);

}

void	Timing_Message(const double cpu1,const double cpu2,
		       const double cpu3,const double cpu4,
		       const double cpustep)
{
   printf("\n");
   printf(">--<>--/~~\\--<>--/~~\\--<>--<\n");
   printf("cpu time in seconds      :( \n");
   printf("prologue     : %le\n", cpu2-cpu1);
   printf("main loop    : %le\n", cpu3-cpu2);
   printf("epilogue     : %le\n", cpu4-cpu3);
   printf("total        : %le\n", cpu4-cpu1);
   printf("cputime/step : %le\n", cpustep);
   printf(">--<>--\\__/--<>--\\__/--<>--<\n");

}

   

void	ELEMENTS_Message(const int nions, char** inames, double** ions)
{
   int i,j;
   int found;
   int nkatms = 0;
   int katms[50];
   int coreq,valenceq;

   printf("  elements involved in the cluster:\n");     

   for (i=0; i<nions; ++i)
   {
      /* Check to make sure that this kind of atom doesnt exist */
      found = 0;
      for (j=0; j<nkatms; ++j)
      {
         if (NameToCharge(inames[i]) == katms[j])
	    found = 1;
      }
      if (!found)
      {
         katms[nkatms] = NameToCharge(inames[i]);
         ++nkatms;
         coreq = NameToCoreCharge(inames[i]);
         valenceq = NameToValenceCharge(inames[i]);

         printf(
	 "      %d: %s      mass no.: %5.1f   core charge: %d   valence charge: %d\n",
	 nkatms, inames[i], NameToMass(inames[i]), coreq,valenceq);
         
      }
   }
   printf("\n");

}
   
   
void	ATOMS_Message(const int nions, char** inames, double** ions)
{
   int i;
   double gcx,gcy,gcz;
   double cmx,cmy,cmz,M;

   printf("  positions of the ions: (bohrs)\n");     

   gcx = 0.0; gcy = 0.0; gcz = 0.0;
   cmx = 0.0; cmy = 0.0; cmz = 0.0;
   M = 0.0;
   for (i=0; i<nions; ++i)
   {
      gcx += ions[i][0];
      gcy += ions[i][1];
      gcz += ions[i][2];
      cmx += ions[i][0]*NameToMass(inames[i]);
      cmy += ions[i][1]*NameToMass(inames[i]);
      cmz += ions[i][2]*NameToMass(inames[i]);
      M   += NameToMass(inames[i]);
      printf(
      "      %d %s\t  (%11.6lf %11.6f %11.6f )\n",
      i+1, inames[i], ions[i][0],ions[i][1],ions[i][2]);
      
   }
   gcx = gcx/nions;
   gcy = gcy/nions;
   gcz = gcz/nions;
   printf("      G.C.\t  (%11.6lf %11.6f %11.6f )\n",gcx,gcy,gcz);
   
   cmx = cmx/M;
   cmy = cmy/M;
   cmz = cmz/M;
   printf(
      "      C.M.\t  (%11.6lf %11.6f %11.6f )\n",cmx,cmy,cmz);
   
   printf("\n");

} /* ATOMS_Message */

void	VELCS_Message(const int nions, char** inames, double** vions)
{
   int i;
   double gcx,gcy,gcz;
   double cmx,cmy,cmz,M;

   printf("  velocity of the ions:  (bohrs)\n");     

   gcx = 0.0; gcy = 0.0; gcz = 0.0;
   cmx = 0.0; cmy = 0.0; cmz = 0.0;
   M = 0.0;
   for (i=0; i<nions; ++i)
   {
      gcx += vions[i][0];
      gcy += vions[i][1];
      gcz += vions[i][2];
      cmx += vions[i][0]*NameToMass(inames[i]);
      cmy += vions[i][1]*NameToMass(inames[i]);
      cmz += vions[i][2]*NameToMass(inames[i]);
      M   += NameToMass(inames[i]);
      printf(
      "      %d %s\t  (%11.6lf %11.6f %11.6f )\n",
      i+1, inames[i], vions[i][0],vions[i][1],vions[i][2]);
      
   }
   gcx = gcx/nions;
   gcy = gcy/nions;
   gcz = gcz/nions;
   printf("      G.C.\t  (%11.6lf %11.6f %11.6f )\n",gcx,gcy,gcz);
   
   cmx = cmx/M;
   cmy = cmy/M;
   cmz = cmz/M;
   printf(
      "      C.M.\t  (%11.6lf %11.6f %11.6f )\n",cmx,cmy,cmz);
   
   printf("\n");

} /* VELCS_Message */






void	ELC_Message(const int nelc, const int ispin, 
		    const int nions, char** inames, double** ions)
{
  int icharge, i;
  printf("  Number of Electrons: %d   ispin: %d\n",nelc,ispin);
  

  icharge = 0;
  for (i=0; i<nions; ++i)
     icharge += NameToValenceCharge(inames[i]);
  icharge = icharge - nelc;
  printf("  Charge of Cluster  : %d\n", icharge);
  
  printf("\n");

}

void	PARAMS_Message(const double dt, const double tole)
{
   printf(
"  technical parameters: time step = %8.3le, tolerence = %8.3le\n\n",
   dt,tole);

   
}

void	PARAMS_TBMD_Message(const double dt, const double scalev)
{
   printf(
"  technical parameters: time step = %8.3le, velocity scaling = %8.3le\n\n",
   dt,scalev);

   
}

