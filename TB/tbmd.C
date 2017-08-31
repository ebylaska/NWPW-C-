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

#define	MAX_IONS	40

#define	Mass_Constant		1822.89
#define Boltzman_Constant	3.16679e-6

double	KE_Ion(const int nions, Ion* vions)
{
   double ke = 0.0;
   for (int i=0; i<nions; ++i)
      ke += 0.5*(vions[i].Mass()*Mass_Constant)* (vions[i].l2norm());

   return ke;
} /* KE_Ion */

double KE_Ion_CM(const int nions, Ion* vions)
{
   double M    = 0.0;
   double cmvx = 0.0;
   double cmvy = 0.0;
   double cmvz = 0.0;
   for (int i=0; i<nions; ++i)
   {
       M    += vions[i].Mass()*Mass_Constant;
       cmvx += vions[i].Mass()*Mass_Constant*vions[i].x();
       cmvy += vions[i].Mass()*Mass_Constant*vions[i].y();
       cmvz += vions[i].Mass()*Mass_Constant*vions[i].z();
   }
   cmvx = cmvx/M;
   cmvy = cmvy/M;
   cmvz = cmvz/M;

   double ke_cm = 0.5*M*(cmvx*cmvx + cmvy*cmvy + cmvz*cmvz);

   return ke_cm;

} /* KE_Ion_CM */
       



main()
{

   int 		nelc, ispin;
   int 		iteration_in, iteration_out;
   int		iterations;

   double	Uin;
   double 	dt,velocity_scale;
   
   double	ke0, ke0_cm;
   double	ke_start, ke_start_cm;
   double	ke,ke_cm;
   double	ke_all, temperature;
   double 	ke_end,ke_end_cm;

   char		tb_parameterization[80];


   int  nions;
   Ion cluster0[MAX_IONS];
   Ion cluster1[MAX_IONS];
   Ion cluster2[MAX_IONS];
   Ion vions[MAX_IONS];
   Ion fions[MAX_IONS];
   double dti[MAX_IONS];
   double h;

   /* used for io */
   double *io_ions[MAX_IONS];
   double *io_vions[MAX_IONS];
   double *io_vions0[MAX_IONS];
   char   *io_names[MAX_IONS];


   double cpu1,cpu2,cpu3,cpu4,cpustep;
   int	  i,j,k;

    seconds(&cpu1);
    Header_TBMD_Message();

    /* Read in the IONIN file */
    {
       ifstream	ion_stream("IONIN");
       Read_xyz(ion_stream,&nions,cluster2);
       ion_stream.close();

       /* set the masses and charges in other ion arrays */
       for (i=0; i<nions; ++i)
       {
          cluster0[i] = cluster2[i];
          cluster1[i] = cluster2[i];
          fions[i]    = cluster2[i];
       }
       /* set the io variables */
       for (i=0; i<nions; ++i)
       {
          io_ions[i]  = cluster2[i].array();
          io_names[i] = cluster2[i].Name();
       }

    }

    /* Read in the VIONIN file */
    {
       ifstream	vion_stream("VIONIN");
       Read_xyz(vion_stream,&nions,cluster1);
       vion_stream.close();

       /* set the io variables */
       for (i=0; i<nions; ++i)
       {
          io_vions[i]   = cluster1[i].array();
       }


    }
  


    /* Read in the CONTROL file */
    {
       ifstream control_stream("CONTROL");
       control_stream >> nelc >> ispin >> Uin;
       control_stream >> dt;
       control_stream >> iteration_in >> iteration_out;
       control_stream >> velocity_scale;
       control_stream >> tb_parameterization;
       control_stream.close();
    }

    /* Read in the tb-parameterization */
    {
      ifstream param_stream(tb_parameterization);
      Read_Interactions(param_stream);
      param_stream.close();
    }


    /* Figure out the kinetic energies and scale velocities */

    ke0    = KE_Ion(nions,cluster1);
    ke0_cm = KE_Ion_CM(nions,cluster1);

    /* scale velocities */
    for (i=0; i<nions; ++i)
       cluster1[i] = velocity_scale*cluster1[i];

    ke_start    = KE_Ion(nions,cluster1);
    ke_start_cm = KE_Ion_CM(nions,cluster1);
   

    /* set up Verlet constants */
    for (i=0; i<nions; ++i)
    {
       double M = cluster2[i].Mass()*Mass_Constant;
       dti[i] = dt*dt/M;
    }
    h = 0.5/dt;
    

    /* output Input Data */

    Input_Data_TBMD_Message(nelc,ispin,dt,velocity_scale,
	nions,io_names,io_ions,io_vions, 
	ke0,ke0_cm,ke_start,ke_start_cm,
	tb_parameterization);


    /* Do the first step */
    /* Initialize the TB routines */

    for (i=0; i<nions; ++i)
    {
       cluster0[i] = cluster1[i];
       cluster1[i] = cluster2[i];
    }
    TB  tb(ispin,nelc,Uin,nions,cluster1,fions);
    tb.Reset();
    tb.Force();

    /* do a newton step */
    for (i=0; i<nions; ++i)
       cluster2[i] = cluster1[i] + dt*cluster0[i] + (0.5*dti[i])*fions[i];


    /* open Motions files */
    ofstream	 motion_stream( "MOTION");
    ofstream	vmotion_stream("VMOTION");

    seconds(&cpu2);
    double a[4];
    int done   = 0;
    iterations = 0;
    Start_Iteration_Message();
    while ((iterations < (iteration_out)) && (!done))
    {
       ++iterations;

       for (int in=1; in<=iteration_in; ++in)
       {
          for (i=0; i<nions; ++i)
          {
             cluster0[i] = cluster1[i];
             cluster1[i] = cluster2[i];
          }
          tb.Reset();
          tb.Force();


          /* do a Verlet step */
          for (i=0; i<nions; ++i)
             cluster2[i] = 2.0*cluster1[i] - cluster0[i] + (dti[i])*fions[i];

	 /* calculate the current velocities */
         for (i=0; i<nions; ++i)
            cluster0[i] = h*(cluster2[i] - cluster0[i]);

         /* calculate the current kinetic energies */
         ke      = KE_Ion(nions,cluster0);
         ke_cm   = KE_Ion_CM(nions,cluster0);
         ke_all += (ke-ke_cm);
       }


       a[0] = tb.Energy() + ke;
       a[1] = tb.Energy();
       a[2] = ke;
       a[3] = temperature;

       Iteration_Message(iterations*iteration_in,4,a);

       /* output to MOTION files */
       Write_xyz(motion_stream,nions,cluster2);
       Write_xyz(vmotion_stream,nions,cluster0);
       

    }
    End_Iteration_Message(done);
    seconds(&cpu3);

    /* close Motions files */
    motion_stream.close();
    vmotion_stream.close();

    /* copy velocity to cluster1 */
    for (i=0; i<nions; ++i)
       cluster1[i] = cluster0[i];


   /* display the output */
   Output_Data_TBMD_Message(tb.Energy(),tb.Core_Energy(),
			tb.Valence_Energy()+tb.U_Energy(),tb.U_Energy(),
			nions,io_names,io_ions,io_vions);

   /* display the eigenvalues */
   Start_Eigenvalue_Message(tb.Spin());
   for (i=0; i<tb.Number_Eigenvalues(); ++i)
      Eigenvalue_Message(i+1,tb.Fill(i), tb(i));

   /* output the IONOUT */
   {
      ofstream	ionout_stream("IONOUT");
      Write_xyz(ionout_stream,nions,cluster2);
      ionout_stream.close();
   }
   /* output the VIONOUT */
   {
      ofstream	vionout_stream("VIONOUT");
      Write_xyz(vionout_stream,nions,cluster1);
      vionout_stream.close();
   }

   /* make the ORBOUT file */
   {
      ofstream	orbout_stream("ORBOUT");
      orbout_stream << tb.Number_Eigenvalues() << " "
         	    << nions << "\n"; 
      orbout_stream << "Pi_Orbs\n";
      for (i=0; i<nions; ++i)
         orbout_stream << cluster2[i].x() << " "
   		       << cluster2[i].y() << " "
   		       << cluster2[i].z() << " 4 \n";
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
   }



     seconds(&cpu4);
     cpustep = (cpu3-cpu2)/(iterations);
     Timing_Message(cpu1,cpu2,cpu3,cpu4,cpustep);

   
    

}
