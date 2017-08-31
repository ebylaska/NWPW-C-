#include	<math.h>
#include	<fstream.h>
#include	<iostream.h>
#include	<string.h>

#include	"seconds.h"
#include	"Message.h"

#include	"Interactions.h"
#include	"Ion.h"
#include	"TB.h"
#include	"xyz_io.h"


#define MAX_IONS	80

using namespace std;

main()
{
   int		i,j,k;
   int 		nelc, ispin;
   int 		iteration_in, iteration_out;
   int		iterations;
   double	Uin;
   double 	dt,tole;
   char		tb_parameterization[80];
   int  nions;
   Ion cluster[MAX_IONS];
   Ion cluster2[MAX_IONS];
   Ion fions[MAX_IONS];

   /* used for io */
   double *io_ions[MAX_IONS];
   char   *io_names[MAX_IONS];

   double cpu1,cpu2,cpu3,cpu4,cpustep;

    seconds(&cpu1);
    Header_Message();

    /* Read in the IONIN file */
    {
       ifstream	ion_stream("IONIN");
       Read_xyz(ion_stream,&nions,cluster);
       ion_stream.close();

       /* set the masses and charges in other ion arrays */
       for (i=0; i<nions; ++i)
       {
          cluster2[i] = cluster[i];
          fions[i]    = cluster[i];
       }
       /* set the io variables */
       for (i=0; i<nions; ++i)
       {
          io_ions[i]  = cluster[i].array();
          io_names[i] = cluster[i].Name();
       }

    }
  


    /* Read in the CONTROL file */
    {
       ifstream control_stream("CONTROL");
       control_stream >> nelc >> ispin >> Uin;
       control_stream >> dt;
       control_stream >> iteration_in >> iteration_out;
       control_stream >> tole;
       control_stream >> tb_parameterization;
       control_stream.close();
    }


    /* Read in the tb-parameterization */
    {
      ifstream param_stream(tb_parameterization);
      Read_Interactions(param_stream);
      param_stream.close();
    }

    Print_Interactions();
    

    /* output Input Data */
    Input_Data_Message(nelc,ispin,dt,tole,
	nions,io_names,io_ions,
	tb_parameterization);

    double a[4];
    double enew,eold;
    iterations = 0;
    seconds(&cpu2);

    /* Initialize the TB routines */
    TB  tb(ispin,nelc,Uin,nions,cluster,fions);
    enew = tb.Energy();


    int done = 0;
    Start_Iteration_Message();
    while ((iterations < (iteration_out*iteration_in)) && (!done))
    {

       for (int in=0; in<iteration_in; ++in)
       {
          tb.Reset();
          tb.Force();
          for (i=0; i<nions; ++i)
          {
             cluster2[i] = cluster[i];
             fions[i]    = fions[i]  * ( dt/sqrt(1822.89*fions[i].Mass()) );
             cluster[i]  = cluster2[i] + fions[i];
          }
       }

       eold = enew;
       enew = tb.Energy();

       a[0] = enew;
       a[1] = tb.Core_Energy();
       a[2] = tb.Valence_Energy() + tb.U_Energy();
       a[3] = fabs(eold - enew);

       iterations += iteration_in;
       Iteration_Message(iterations,4,a);

       if (a[3] < tole) done = 1;
    }
    End_Iteration_Message(done);
    seconds(&cpu3);

   /* display the output */
   Output_Data_Message(tb.Energy(),tb.Core_Energy(),
			tb.Valence_Energy()+tb.U_Energy(),tb.U_Energy(),
			nions,io_names,io_ions);

   /* display the eigenvalues */
   Start_Eigenvalue_Message(tb.Spin());
   for (i=0; i<tb.Number_Eigenvalues(); ++i)
      Eigenvalue_Message(i+1,tb.Fill(i), tb(i));

   /* output the IONOUT */
   ofstream	ionout_stream("IONOUT");
   Write_xyz(ionout_stream,nions,cluster);
   ionout_stream.close();

   /* make the ORBOUT file */
   ofstream	orbout_stream("ORBOUT");
   orbout_stream << tb.Number_Eigenvalues() << " "
		 << nions << "\n"; 
   orbout_stream << "Pi_Orbs\n";
   for (i=0; i<nions; ++i)
      orbout_stream << cluster[i].x() << " "
		    << cluster[i].y() << " "
		    << cluster[i].z() << " 4 \n";
   for (i=0; i<tb.Number_Eigenvalues(); ++i)
   {
      for (j=0; j<nions; ++j)
      {
          for (k=0; k<4; ++k)
             orbout_stream << tb(k+j*4,i) << " ";
         orbout_stream <<"\n";
      }
   }
   orbout_stream.close();



     seconds(&cpu4);
     cpustep = (cpu3-cpu2)/(iterations);
     Timing_Message(cpu1,cpu2,cpu3,cpu4,cpustep);

   
    

}
