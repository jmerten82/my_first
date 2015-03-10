#include <iostream>
#include <ctime>
#include <cmath>
#include <gsl/gsl_matrix.h>
#include "mpi.h"
#include "options.h"
#include "rw_fits.h"
#include "mpi_comm.h"
#include "mpi_coreroutines.h"
#include "mpi_datacontainer.h"
#include "interpol.h"
#include "soph_math.h"


using namespace std;

int main(int argc, char *argv[])
{

  /*
    Starting up and some definitions.
  */

  int my_rank; //process rank
  int p; //total numver of processes


  //init MPI and read process info

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&p);



  //doing some timing and shiny startup screen, since we only want to see that once
  //we do it on the head process.

  long starttime;
  long stoptime;

  if(my_rank == 0)
    {
      cout <<"#------------------SaWLens 1.2 Reconstruction tool------------------#" <<endl;
      cout <<endl;
      cout <<"Running on: " <<p <<" processes." <<endl;
      time(&starttime);
    }

  //Read the configuration file on all processes.

  ReconstructionOptions options(argv[1]);

  //Defining some auxilliary quantities which are needed during 
  //reconstruction.

  int x_dim; //Pixel dimensions 
  int y_dim;
  int x_dimplusone;
  int y_dimplusone;
  gsl_matrix_int *mask;  //field mask file, more or less defining the reconstr.
  gsl_matrix *interpolconv; //will become the interpolated convergence map
  gsl_matrix *interpolshear1; //etc.
  gsl_matrix *interpolshear2;
  gsl_matrix *interpolf1;
  gsl_matrix *interpolf2;
  gsl_matrix *interpolg1;
  gsl_matrix *interpolg2;
  gsl_vector *control; //will check for the convergence change during rec. 
  gsl_vector *pot; //Carries some intermediate results through the inner-level
                  //iterations

  int ilevel; //iteration index needed later




  /**
    Start of the two-level reconstruction process.
  **/



 
  //"outer-level" iteration loop, see arXix:0806.1967
  for(int olevel = 0; olevel < options.shownumsteps(); olevel++)
    {
      MPI_Barrier(MPI_COMM_WORLD);

      //Getting grid sizes
      if(options.showbool("shear"))
	{
	  x_dim = read_intheader(options.show_filename("shearinput",olevel),"NAXIS1");
	  y_dim = read_intheader(options.show_filename("shearinput",olevel),"NAXIS2");
	  if(olevel != options.shownumsteps()-1)
	    {
	      x_dimplusone = read_intheader(options.show_filename("shearinput",olevel+1),"NAXIS1");
	      y_dimplusone = read_intheader(options.show_filename("shearinput",olevel+1),"NAXIS2");
	    }
	}
      else
	{
	  x_dim = read_intheader(options.show_filename("flexioninput",olevel),"NAXIS1");
	  y_dim = read_intheader(options.show_filename("flexioninput",olevel),"NAXIS2");
	  if(olevel != options.shownumsteps()-1)
	    {
	      x_dimplusone = read_intheader(options.show_filename("flexioninput",olevel+1),"NAXIS1");
	      y_dimplusone = read_intheader(options.show_filename("flexioninput",olevel+1),"NAXIS2");
	    }
	}

      if(my_rank == 0)
	{
	  cout <<endl;
	  cout <<"Outer-level iteration started." <<endl;
	  cout <<"Reconstruction resolution: " <<x_dim <<"x" <<y_dim <<endl;
	}

      //Allocating, yes C++ is RowMajor which is not really helpful while 
      //thinking about map-like grids

      mask = gsl_matrix_int_calloc(y_dim,x_dim);
      

      //Reading the mask for this resolution and defining the grid,
      //disk access only on leading process and then sending
      if(my_rank == 0)
	{
	  if(options.showbool("shear"))
	    {
	      read_imgeint(options.show_filename("shearinput",olevel),"field_mask",mask);
	    }
	  else
	    {
	      read_imgeint(options.show_filename("flexioninput",olevel),"field_mask",mask);
	    }
	}
      send_gsl_toworld(mask);
      FinDifGrid grid(mask,1.0,true);
      //allocating helper quantities, actually pot will often carry the result
      control = gsl_vector_calloc(grid.showint("fieldpixels"));
      pot = gsl_vector_calloc(grid.showint("fieldpixels"));




      //Reading the data for the reconstruction and sending to all processes

      MPIDataContainer container(options,grid,my_rank,p);
      container.readData(options,grid,olevel);
      //sync before MPI method
      MPI_Barrier(MPI_COMM_WORLD);
      container.sendData();

      //Building a reference for regularisation

      GalaxyCluster reference(mask,1.0);
      //at the first step we don't have interpolated results
      if(olevel != 0)
	{
	  reference.readmap(interpolconv,"convergence");
	  reference.readmap(interpolshear1,"shear1");
	  reference.readmap(interpolshear2,"shear2");
	  reference.readmap(interpolf1,"f1");
	  reference.readmap(interpolf2,"f2");
	  reference.readmap(interpolg1,"g1");
	  reference.readmap(interpolg2,"g2");
	  reference.buildwithoutpot();
	}
      //else all values are set to 0 which is also true for convergence
      //here we effectively set the often discussed flat convergence prior


      //clearing memory if we allocated it
      if(olevel != 0)
	{
	  gsl_matrix_free(interpolconv);
	  gsl_matrix_free(interpolshear1);
	  gsl_matrix_free(interpolshear2);
	  gsl_matrix_free(interpolf1);
	  gsl_matrix_free(interpolf2);
	  gsl_matrix_free(interpolg1);
	  gsl_matrix_free(interpolg2);
	}

    
      //"inner-level" iteration loop
      double threshold = 100.0;
      ilevel = 0;
      while((ilevel < options.showint("iterations")) && (threshold > options.showdouble("threshold")))
	{

	  MPI_Barrier(MPI_COMM_WORLD);
	  if(my_rank == 0)
	    {
	      cout <<endl;
	      cout <<"Inner level iteration: " <<ilevel+1 <<endl;
	      cout <<endl;
	    }
	  MPI_Barrier(MPI_COMM_WORLD);
	  //Building a cluster which contains the reconstruction values during
	  //the inner level iterations
	  GalaxyCluster cluster(mask,1.0);
	  if(ilevel != 0)
	    {
	      cluster.readgridvector(pot,"pot");
	      cluster.buildfrompot();
	    }
	  else
	    {
	      cluster.readgridvector(reference.data("convergence"),"convergence");
	      cluster.readgridvector(reference.data("shear1"),"shear1");
	      cluster.readgridvector(reference.data("shear2"),"shear2");
	      cluster.readgridvector(reference.data("f1"),"f1");
	      cluster.readgridvector(reference.data("f2"),"f2");
	      cluster.readgridvector(reference.data("g1"),"g1");
	      cluster.readgridvector(reference.data("g2"),"g2");
	      cluster.readgridvector(reference.data("jacdet"),"jacdet");
	    }
	  //keeping the former convergence value
	  gsl_vector_memcpy(control,cluster.data("convergence"));
	  //setting reconstruction quantities to zero on all processes
	  container.resetResult();
	  //Building prefactors
	  if(!options.showbool("weaklimit"))
	      {
		  MPIWeakFactor(options,cluster,container);
	      }
	  MPIStrongFactor(options,cluster,container);
	  MPI_Blk(options,reference,container,p,my_rank);
	  //this sync is actually just for output style reasons, not REALLY necessary
	  MPI_Barrier(MPI_COMM_WORLD);
	  MPI_Vl(options,reference,container,p,my_rank);
	  //sync before talking
	  MPI_Barrier(MPI_COMM_WORLD);
	  container.sendResult();
	  MPI_Barrier(MPI_COMM_WORLD);
	  //solving linear system
	  if(my_rank == 0)
	    {
	      solve_gsl(container.showdatamatrix("finalcoeff"),container.showdatavector("finaldatavector"));
	      gsl_vector_memcpy(pot,container.showdatavector("finaldatavector"));
	    }
	  //sync before talking
	  MPI_Barrier(MPI_COMM_WORLD);
	  send_gsl_toworld(pot);


	  cluster.readgridvector(pot,"pot");
	  cluster.buildfrompot();
	  if(my_rank == 0)
	    {
	      cout <<endl;
	    }
	  threshold = change(cluster.data("convergence"),control);
	  if(my_rank == 0)
	    {
	      cout <<"Maximum change in convergence: " <<threshold <<endl;
	    }
	  if(options.showbool("weaklimit"))
	      {
		  ilevel = options.showint("iterations"); //just to make sure
	      }
	  ilevel++;
	}
    

      

      //Writing the result to FITS
      if(my_rank == 0)
      {
	GalaxyCluster output(mask,1.0);
	output.readgridvector(pot,"pot");
	output.buildfrompot();
	output.writetofits(options,olevel);
      }

      //Performing interpolations, including map allocations
      if(olevel != options.shownumsteps()-1)
	{
	  interpolconv = gsl_matrix_calloc(y_dimplusone,x_dimplusone);
	  interpolshear1 = gsl_matrix_calloc(y_dimplusone,x_dimplusone);
	  interpolshear2 = gsl_matrix_calloc(y_dimplusone,x_dimplusone);
	  interpolf1 = gsl_matrix_calloc(y_dimplusone,x_dimplusone);
	  interpolf2 = gsl_matrix_calloc(y_dimplusone,x_dimplusone);
	  interpolg1 = gsl_matrix_calloc(y_dimplusone,x_dimplusone);
	  interpolg2 = gsl_matrix_calloc(y_dimplusone,x_dimplusone);

	  if(my_rank == 0)
	    {
	      GalaxyCluster interpol(mask,1.0);
	      interpol.readgridvector(pot,"pot");
	      interpol.buildfrompot();
	      gsl_matrix * nofuture1 = gsl_matrix_calloc(y_dim,x_dim);
	      interpol.writemap(nofuture1,"convergence");
	      gsl_matrix * nofuture2 = gsl_matrix_calloc(y_dim,x_dim);
	      interpol.writemap(nofuture2,"shear1");
	      gsl_matrix * nofuture3 = gsl_matrix_calloc(y_dim,x_dim);
	      interpol.writemap(nofuture3,"shear2");
	      gsl_matrix * nofuture4 = gsl_matrix_calloc(y_dim,x_dim);
	      interpol.writemap(nofuture4,"f1");
	      gsl_matrix * nofuture5 = gsl_matrix_calloc(y_dim,x_dim);
	      interpol.writemap(nofuture5,"f2");
	      gsl_matrix * nofuture6 = gsl_matrix_calloc(y_dim,x_dim);
	      interpol.writemap(nofuture6,"g1");
	      gsl_matrix * nofuture7 = gsl_matrix_calloc(y_dim,x_dim);
	      interpol.writemap(nofuture7,"g2");
	      cspline_morecleverinterpol(nofuture1,mask,interpolconv);
	      cspline_morecleverinterpol(nofuture2,mask,interpolshear1);
	      cspline_morecleverinterpol(nofuture3,mask,interpolshear2);
	      cspline_morecleverinterpol(nofuture4,mask,interpolf1);
	      cspline_morecleverinterpol(nofuture5,mask,interpolf2);
	      cspline_morecleverinterpol(nofuture6,mask,interpolg1);
	      cspline_morecleverinterpol(nofuture7,mask,interpolg2);

	      gsl_matrix_free(nofuture1);
	      gsl_matrix_free(nofuture2);
	      gsl_matrix_free(nofuture3);
	      gsl_matrix_free(nofuture4);
	      gsl_matrix_free(nofuture5);
	      gsl_matrix_free(nofuture6);
	      gsl_matrix_free(nofuture7);
	    }
	  //sync before talking
	  MPI_Barrier(MPI_COMM_WORLD);
	  send_gsl_toworld(interpolconv);
	  send_gsl_toworld(interpolshear1);
	  send_gsl_toworld(interpolshear2);
	  send_gsl_toworld(interpolf1);
	  send_gsl_toworld(interpolf2);
	  send_gsl_toworld(interpolf1);
	  send_gsl_toworld(interpolf2);
	  send_gsl_toworld(interpolg1);
	  send_gsl_toworld(interpolg2);
	  MPI_Barrier(MPI_COMM_WORLD);
	}
      //clearing aux quantities
      gsl_matrix_int_free(mask);
      gsl_vector_free(pot);
      gsl_vector_free(control);
    }

  if(my_rank == 0)
    {
      time(&stoptime);
      cout <<endl;
      cout <<"Total runtime: "<<stoptime-starttime <<" sec." <<endl;
    }


  MPI_Finalize();

  return 0;

}

            

		      
