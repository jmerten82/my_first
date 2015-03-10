#include <iostream>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include "rw_fits.h"
#include "mpi.h"
#include "mpi_coreroutines.h"
#include "mpi_comm.h"
#include "galaxycluster.h"
#include "reconstruction_kernel.h"
#include "nrinterpol.h"
#include "util.h"
#include "options.h"

using namespace std;

int main(int argc, char *argv[])
{

  int my_rank;
  int p;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&p);

  long starttime;
  long stoptime;
  string infile;
  string outfile;

  BootstrapReconstructionOptions options1(argv[1]);

  for(int b = options1.showbootstrapindex(0); b <= options1.showbootstrapindex(1); b++)
    {

      gsl_vector *ellip1 = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_vector *ellip2 = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_matrix *Cij1 = gsl_matrix_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0),options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_matrix *Cij2 = gsl_matrix_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0),options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_matrix *Fij1 = gsl_matrix_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0),options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_matrix *Fij2 = gsl_matrix_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0),options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_vector *wredshift = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_vector *recpot = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_vector *recconv = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_vector *recshear1 = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_vector *recshear2 = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_vector *refconv = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_vector *refshear1 = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_vector *refshear2 = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_vector *eta = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_matrix *proccoeff = gsl_matrix_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0),options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_matrix *calccoeff = gsl_matrix_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0),options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_vector *procdata = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_vector *calcdata = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      double strongsigma;
      double control;
      gsl_vector *nofuture1 = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_vector *control1 = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_vector *control2 = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_vector *dummy = gsl_vector_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));
      gsl_vector_int *dummy2 = gsl_vector_int_calloc(options1.showstepbyindex(0)*options1.showstepbyindex(0));


      ostringstream in;
      ostringstream out;

      in <<options1.showfilename(0) <<"/ellip_realisation" <<b <<"_" <<options1.showstepbyindex(0) <<".fits";
      out <<options1.showfilename(1) <<"/rec_realisation" <<b <<"_" <<options1.showstepbyindex(0) <<".fits";

      infile = in.str();
      outfile = out.str();

      if(my_rank == 0)
	{

	  cout <<"------------SAW bootstrapping reconstruction routine------------------" <<endl;
	  cout <<"Processing realisation: " <<b <<endl;

	  read_pimg(infile,ellip1);
	  read_imge(infile,"av_ellip2",ellip2);
	  read_imge(infile,"ellip_covariance1",Cij1);
	  read_imge(infile,"ellip_covariance2",Cij2);
	  gsl_vector_set_all(eta,options1.showreg());
	  gsl_vector_set_all(wredshift,cosmicweight(options1.showredshift(1),options1.showredshift(0)));
	}

      send_gsl_toworld_new(ellip1);
      send_gsl_toworld_new(ellip2);
      send_gsl_toworld_new(Cij1);
      send_gsl_toworld_new(Cij2);
      send_gsl_toworld_new(eta);
      send_gsl_toworld_new(wredshift);
      MPI_Barrier(MPI_Comm MPI_COMM_WORLD);

      if(my_rank == 0)
	{
	  time(&starttime);
	  GalaxyCluster cluster1(options1.showresolution(0),options1.showresolution(0),1.0);
	  weak_init(options1.showstepbyindex(0),options1.showstepbyindex(0),1.0,dummy2,ellip1,ellip2,Cij1,Cij2,refconv,refshear1,refshear2,eta,dummy,wredshift,dummy,recpot);
	  cluster1.readfrom("pot",recpot);
	  cluster1.buildfrompot();
	}

      gsl_vector_memcpy(nofuture1,recpot);
      MPI_Barrier(MPI_Comm MPI_COMM_WORLD);

      if(my_rank == 0)
	{
	  cout <<"Starting Outer-Level iteration." <<endl;
	}
      
      for(int i = 0; i < options1.shownumofsteps(); i++)
	{
	  MPI_Barrier(MPI_Comm MPI_COMM_WORLD);
	  if(my_rank == 0)
	    {
	      cout <<"Reconstruction resolution: " <<options1.showstepbyindex(i) <<"x" <<options1.showstepbyindex(i) <<endl;
	    }


	      ostringstream in;
	      ostringstream out;
	      in <<options1.showfilename(0) <<"/ellip_realisation" <<b <<"_" <<options1.showstepbyindex(i) <<".fits";
	      out <<options1.showfilename(1) <<"/rec_realisation" <<b <<"_" <<options1.showstepbyindex(i) <<".fits";
	      infile = in.str();
	      outfile = out.str();

	      ellip1 = gsl_vector_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      ellip2 = gsl_vector_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      dummy = gsl_vector_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      dummy2 = gsl_vector_int_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      Cij1 = gsl_matrix_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i),options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      Cij2 = gsl_matrix_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i),options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      eta = gsl_vector_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      wredshift = gsl_vector_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      recpot = gsl_vector_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      recconv = gsl_vector_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      recshear1 = gsl_vector_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      recshear2 = gsl_vector_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      refconv = gsl_vector_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      refshear1 = gsl_vector_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      refshear2 = gsl_vector_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      control1 = gsl_vector_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      control2 = gsl_vector_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      Fij1 = gsl_matrix_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i),options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      Fij2 = gsl_matrix_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i),options1.showstepbyindex(i)*options1.showstepbyindex(i));
	      MPI_Barrier(MPI_Comm MPI_COMM_WORLD);

	      if(my_rank == 0)
		{
		  gsl_vector_memcpy(recpot,nofuture1);
		  GalaxyCluster cluster1(options1.showstepbyindex(i),options1.showstepbyindex(i),1.0);
		  cluster1.readfrom("pot",recpot);
		  cluster1.buildfrompot();
		  //cluster1.masssheetnormalise();
		  cluster1.writeto("conv",recconv);
		  cluster1.writeto("shear1",recshear1);
		  cluster1.writeto("shear2",recshear2);
		  gsl_vector_free(nofuture1);
		}
	      MPI_Barrier(MPI_Comm MPI_COMM_WORLD);
	      
	      send_gsl_toworld_new(recconv);
	      gsl_vector_memcpy(refconv,recconv);
	      send_gsl_toworld_new(recshear1);
	      gsl_vector_memcpy(refshear1,recshear1);
	      send_gsl_toworld_new(recshear2);
	      gsl_vector_memcpy(refshear2,recshear2);
	      MPI_Barrier(MPI_Comm MPI_COMM_WORLD);

	      if(my_rank == 0)
		{
		  read_pimg(infile,ellip1);
		  read_imge(infile,"av_ellip2",ellip2);
		  read_imge(infile,"ellip_covariance1",Cij1);
		  read_imge(infile,"ellip_covariance2",Cij2);
		  gsl_vector_set_all(wredshift,cosmicweight(options1.showredshift(1),options1.showredshift(0)));
		  gsl_vector_set_all(eta,options1.showreg());
		}


	      send_gsl_toworld_new(ellip1);
	      send_gsl_toworld_new(ellip2);
	      send_gsl_toworld_new(Cij1);
	      send_gsl_toworld_new(Cij2);
	      send_gsl_toworld_new(wredshift);
	      send_gsl_toworld_new(eta);
	      control = 1.0;
	      MPI_Barrier(MPI_Comm MPI_COMM_WORLD);

	      if(my_rank == 0)
		{
		  cout <<"Starting Inner-Level iteration." <<endl;
		}
	      
	      for(int j = 0; j < 3; j++)
		{
		  if(my_rank == 0)
		    {
		      cout <<"Iteration step: " <<j+1 <<endl;
		    }
		  
		  gsl_vector_memcpy(control1,recconv);
		  Weakfactor(options1.showstepbyindex(i)*options1.showstepbyindex(i),options1.showstepbyindex(i)*options1.showstepbyindex(i),recconv,wredshift,Fij1,Fij2,Cij1,Cij2);
		  MPI_Barrier(MPI_Comm MPI_COMM_WORLD);


		  procdata = gsl_vector_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
		  calcdata = gsl_vector_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i));
		  proccoeff = gsl_matrix_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i),options1.showstepbyindex(i)*options1.showstepbyindex(i));
		  calccoeff = gsl_matrix_calloc(options1.showstepbyindex(i)*options1.showstepbyindex(i),options1.showstepbyindex(i)*options1.showstepbyindex(i));
		  Shear_data_vector(options1.showstepbyindex(i),options1.showstepbyindex(i),1.0,dummy2,ellip1,ellip2,Fij1,Fij2,refconv,refshear1,refshear2,eta,dummy,wredshift,dummy,procdata,my_rank);
		  MPI_Barrier(MPI_Comm MPI_COMM_WORLD);
		  Shear_coeff_matrix(options1.showstepbyindex(i),options1.showstepbyindex(i),1.0,dummy2,ellip1,ellip2,Fij1,Fij2,eta,dummy,wredshift,dummy,proccoeff,my_rank);
		  MPI_Barrier(MPI_Comm MPI_COMM_WORLD);
		  
		  recv_gsl_fromworld_new(procdata,calcdata);
		  MPI_Barrier(MPI_Comm MPI_COMM_WORLD);
		  recv_gsl_fromworld_new(proccoeff,calccoeff);
		  MPI_Barrier(MPI_Comm MPI_COMM_WORLD);

		  if(my_rank == 0)
		    {
		      Solve_system(calccoeff,calcdata,recpot,options1.showstepbyindex(i)*options1.showstepbyindex(i));
		      GalaxyCluster cluster1(options1.showstepbyindex(i),options1.showstepbyindex(i),1.0);
		      cluster1.readfrom("pot",recpot);
		      cluster1.buildfrompot();
		      //cluster1.masssheetnormalise();
		      cluster1.writeto("conv",recconv);
		      cluster1.writeto("shear1",recshear1);
		      cluster1.writeto("shear2",recshear2);
		      gsl_vector_memcpy(control2,recconv);
		      control = change(control1,control2,cluster1.showdim());
		      cout <<"Maximum change in convergence: " <<control <<endl;
		    }
		  
		  send_gsl_toworld_new(recconv);
		  send_gsl_toworld_new(recshear1);
		  send_gsl_toworld_new(recshear2);
		  send_double_toworld(control,p,my_rank);
		  MPI_Barrier(MPI_Comm MPI_COMM_WORLD);

		  gsl_vector_free(procdata);
		  gsl_vector_free(calcdata);
		  gsl_matrix_free(proccoeff);
		  gsl_matrix_free(calccoeff);
		  
		}
	      
	      if(my_rank == 0)
		{
		  GalaxyCluster cluster1(options1.showstepbyindex(i),options1.showstepbyindex(i),1.0);
		  cluster1.readfrom("pot",recpot);
		  cluster1.buildfrompot();
		  //cluster1.masssheetnormalise();
		  cluster1.writeto("conv",recconv);
		  cluster1.writeto("shear1",recshear1);
		  cluster1.writeto("shear2",recshear2);
		  if(i == options1.shownumofsteps()-1)
		    {
		      cluster1.writetofits("all",outfile);
		    }

		  
		}
	      MPI_Barrier(MPI_Comm MPI_COMM_WORLD);

	      if(i != (options1.shownumofsteps()-1))
		{
		  
		  nofuture1 = gsl_vector_calloc(options1.showstepbyindex(i+1)*options1.showstepbyindex(i+1));
		  if(my_rank == 0)
		    {
		      smoothinterpol(recpot,options1.showstepbyindex(i),options1.showstepbyindex(i),options1.showstepbyindex(i+1)-options1.showstepbyindex(i),nofuture1);
		    }
		}
	      
	      MPI_Barrier(MPI_Comm MPI_COMM_WORLD);

	      gsl_vector_free(ellip1);
	      gsl_vector_free(ellip2);
	      gsl_vector_free(wredshift);
	      gsl_vector_free(recpot);
	      gsl_vector_free(recconv);
	      gsl_vector_free(recshear1);
	      gsl_vector_free(recshear2);
	      gsl_vector_free(refconv);
	      gsl_vector_free(refshear1);
	      gsl_vector_free(refshear2);
	      gsl_vector_free(control1);
	      gsl_vector_free(control2);
	      gsl_vector_free(eta);
	      //gsl_vector_free(nofuture1);
	      gsl_vector_free(dummy);
	      gsl_vector_int_free(dummy2);
	      
	      gsl_matrix_free(Cij1);
	      gsl_matrix_free(Cij2);
	      gsl_matrix_free(Fij1);
	      gsl_matrix_free(Fij2);
	      
	}

    }

  MPI_Finalize();
  
      
  
  return 0;
      
}



	  
 


      
