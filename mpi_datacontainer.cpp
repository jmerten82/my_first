#include "mpi_datacontainer.h"

using namespace std;

MPIDataContainer::MPIDataContainer(ReconstructionOptions &options,FinDifGrid &grid1,int my_rankgiven,int pgiven)
{
  my_rank = my_rankgiven;
  p = pgiven;
  
  MPI_Array_Manager array1(my_rank,p,grid1.showint("fieldpixels"),0);
  shear = options.showbool("shear");
  flexion = options.showbool("flexion");
  ccurveflag = options.showbool("ccurve");
  msystems = options.showbool("msystems");
  created = array1.greenlight();
  root = array1.showroot();

  my_rank = my_rankgiven;
  p = pgiven;
  dim = grid1.showint("fieldpixels");

  if(array1.greenlight())
    {
      coeff = gsl_matrix_calloc(array1.showsize(),grid1.showint("fieldpixels"));
      datavector = gsl_vector_calloc(array1.showsize());
    }

  if(root)
    {
      finalcoeff = gsl_matrix_calloc(grid1.showint("fieldpixels"),grid1.showint("fieldpixels"));
      finaldatavector = gsl_vector_calloc(grid1.showint("fieldpixels"));
    }
  else
    {
      finalcoeff = gsl_matrix_calloc(1,1);
      finaldatavector = gsl_vector_calloc(1);
    }

  if(shear)
    {
      ellip1 = gsl_vector_calloc(grid1.showint("fieldpixels"));
      ellip2 = gsl_vector_calloc(grid1.showint("fieldpixels"));
      ellip1covariance = gsl_matrix_calloc(grid1.showint("fieldpixels"),grid1.showint("fieldpixels"));
      ellip2covariance = gsl_matrix_calloc(grid1.showint("fieldpixels"),grid1.showint("fieldpixels"));
      Fij1 = gsl_matrix_calloc(grid1.showint("fieldpixels"),grid1.showint("fieldpixels"));
      Fij2 = gsl_matrix_calloc(grid1.showint("fieldpixels"),grid1.showint("fieldpixels"));
      shear_redshift = gsl_vector_calloc(grid1.showint("fieldpixels"));
    }
  if(flexion)
    {
      f1 = gsl_vector_calloc(grid1.showint("fieldpixels"));
      f2 = gsl_vector_calloc(grid1.showint("fieldpixels"));
      f1covariance = gsl_matrix_calloc(grid1.showint("fieldpixels"),grid1.showint("fieldpixels"));
      f2covariance = gsl_matrix_calloc(grid1.showint("fieldpixels"),grid1.showint("fieldpixels"));
      g1 = gsl_vector_calloc(grid1.showint("fieldpixels"));
      g2 = gsl_vector_calloc(grid1.showint("fieldpixels"));
      g1covariance = gsl_matrix_calloc(grid1.showint("fieldpixels"),grid1.showint("fieldpixels"));
      g2covariance = gsl_matrix_calloc(grid1.showint("fieldpixels"),grid1.showint("fieldpixels"));
      flexion_redshift = gsl_vector_calloc(grid1.showint("fieldpixels"));

    }
  if(ccurveflag)
    {
      ccurve = gsl_vector_int_calloc(grid1.showint("fieldpixels"));
      ccurve_redshift = gsl_vector_calloc(grid1.showint("fieldpixels"));
      ccurveerror = gsl_vector_calloc(grid1.showint("fieldpixels"));
      strongfactorBlk = gsl_vector_calloc(grid1.showint("fieldpixels"));
      strongfactorVl = gsl_vector_calloc(grid1.showint("fieldpixels"));
      
    }
  
  if(msystems)
    {
      msysteminfo = gsl_matrix_calloc(18,options.showint("nummsystems"));
    }
  
}


MPIDataContainer::~MPIDataContainer()
{
  if(created)
    {      
      gsl_matrix_free(coeff);
      gsl_vector_free(datavector);
    }
  gsl_matrix_free(finalcoeff);
  gsl_vector_free(finaldatavector);
  if(shear)
    {
      gsl_vector_free(ellip1);
      gsl_vector_free(ellip2);
      gsl_matrix_free(ellip1covariance);
      gsl_matrix_free(ellip2covariance);
      gsl_matrix_free(Fij1);
      gsl_matrix_free(Fij2);
      gsl_vector_free(shear_redshift);
    }
  if(flexion)
    {
      gsl_vector_free(f1);
      gsl_vector_free(f2);
      gsl_matrix_free(f1covariance);
      gsl_matrix_free(f2covariance);
      gsl_vector_free(g1);
      gsl_vector_free(g2);
      gsl_matrix_free(g1covariance);
      gsl_matrix_free(g2covariance);
      gsl_vector_free(flexion_redshift);
    }
  if(ccurveflag)
    {
      gsl_vector_int_free(ccurve);
      gsl_vector_free(ccurve_redshift);
      gsl_vector_free(ccurveerror);
      gsl_vector_free(strongfactorBlk);
      gsl_vector_free(strongfactorVl);

    }
  if(msystems)
    {
      gsl_matrix_free(msysteminfo);
    }
}

void MPIDataContainer::readData(ReconstructionOptions &options,FinDifGrid &grid, int iterationindex)
{

  if(my_rank == 0)
    {
      gsl_matrix *dummymap = gsl_matrix_calloc(grid.showint("y_dim"),grid.showint("x_dim"));
      gsl_matrix_int *dummymapint = gsl_matrix_int_calloc(grid.showint("y_dim"),grid.showint("x_dim"));

      if(shear)
	{
	  read_pimg(options.show_filename("shearinput",iterationindex),dummymap);
	  grid.maptogridvector(dummymap,ellip1);
	  read_imge(options.show_filename("shearinput",iterationindex),"mean_ellip2",dummymap);
	  grid.maptogridvector(dummymap,ellip2);
	  read_imge(options.show_filename("shearinput",iterationindex),"ellip1_covariance",ellip1covariance);
	  read_imge(options.show_filename("shearinput",iterationindex),"ellip2_covariance",ellip2covariance);
	  gsl_vector_set_all(shear_redshift,options.showdouble("shearz")); 
	}
      
      if(flexion)
	{
	  read_pimg(options.show_filename("flexioninput",iterationindex),dummymap);
	  grid.maptogridvector(dummymap,f1);
	  read_imge(options.show_filename("flexioninput",iterationindex),"mean_f2",dummymap);
	  grid.maptogridvector(dummymap,f2);
	  read_imge(options.show_filename("flexioninput",iterationindex),"mean_g1",dummymap);
	  grid.maptogridvector(dummymap,g1);
	  read_imge(options.show_filename("flexioninput",iterationindex),"mean_g2",dummymap);
	  grid.maptogridvector(dummymap,g2);
	  read_imge(options.show_filename("flexioninput",iterationindex),"f1_covariance",f1covariance);
	  read_imge(options.show_filename("flexioninput",iterationindex),"f2_covariance",f2covariance);
	  read_imge(options.show_filename("flexioninput",iterationindex),"g1_covariance",g1covariance);
	  read_imge(options.show_filename("flexioninput",iterationindex),"g2_covariance",g2covariance);
	  gsl_vector_set_all(flexion_redshift,options.showdouble("flexionz"));
	}
	  
      if(ccurveflag)
	{
	  read_pimgint(options.show_filename("ccurveinput",iterationindex),dummymapint);
	  grid.maptogridvector(dummymapint,ccurve);
	  read_imge(options.show_filename("ccurveinput",iterationindex),"redshift_info",dummymap);
	  grid.maptogridvector(dummymap,ccurve_redshift);
	  read_imge(options.show_filename("ccurveinput",iterationindex),"error",dummymap);
	  grid.maptogridvector(dummymap,ccurveerror);
	} 

      if(msystems)
	{
	  ReadMsystemInfo(options.show_filename("msystemsinput",iterationindex),msysteminfo);
	}
      gsl_matrix_free(dummymap);
      gsl_matrix_int_free(dummymapint);

    }
}

void MPIDataContainer::sendData()
{
  MPI_Barrier(MPI_COMM_WORLD);
  if(shear)
    {
      send_gsl_toworld(ellip1);
      send_gsl_toworld(ellip2);
      send_gsl_toworld(ellip1covariance);
      send_gsl_toworld(ellip2covariance);
      send_gsl_toworld(Fij1);
      send_gsl_toworld(Fij2);
      send_gsl_toworld(shear_redshift);
    }
  if(flexion)
    {
      send_gsl_toworld(f1);
      send_gsl_toworld(f2);
      send_gsl_toworld(f1covariance);
      send_gsl_toworld(f2covariance);
      send_gsl_toworld(g1);
      send_gsl_toworld(g2);
      send_gsl_toworld(g1covariance);
      send_gsl_toworld(g2covariance);
      send_gsl_toworld(flexion_redshift);
    }
  if(ccurveflag)
    {
      send_gsl_toworld(ccurve);
      send_gsl_toworld(ccurve_redshift);
      send_gsl_toworld(ccurveerror);

    }
  if(msystems)
    {
      send_gsl_toworld(msysteminfo);
    }
}

void MPIDataContainer::sendResult()
{
  MPI_Array_Manager array1(my_rank,p,dim,0);
  array1.createcomm();
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(array1.greenlight())
    {
      recv_gsl_fromworld(coeff,finalcoeff,array1);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  if(array1.greenlight())
    {
  recv_gsl_fromworld(datavector,finaldatavector,array1);
    }
  MPI_Barrier(MPI_COMM_WORLD);
   

}
void MPIDataContainer::resetResult()
{
  if(shear)
    {
      gsl_matrix_set_all(Fij1,0.0);
      gsl_matrix_set_all(Fij2,0.0);
    }
  if(ccurveflag)
    {
      gsl_vector_set_all(strongfactorBlk,0.0);
      gsl_vector_set_all(strongfactorVl,0.0);
    }
  gsl_matrix_set_all(finalcoeff,0.0);
  gsl_vector_set_all(finaldatavector,0.0);
  if(created)
    {
      gsl_matrix_set_all(coeff,0.0);
      gsl_vector_set_all(datavector,0.0);
    }
}

gsl_vector* MPIDataContainer::showdatavector(string selection)
{
  if(selection == "ellip1")
    {
      return ellip1;
    }
  else if(selection =="ellip2")
    {
      return ellip2;
    }
  else if(selection == "shear_redshift")
    {
      return shear_redshift;
    }
  else if(selection == "f1")
    {
      return f1;
    }
  else if(selection == "f2")
    {
      return f2;
    }
  else if(selection == "g1")
    {
      return g1;
    }
  else if(selection == "g2")
    {
      return g2;
    }
  else if(selection == "flexion_redshift")
    {
      return flexion_redshift;
    }
  else if(selection == "ccurve_redshift")
    {
      return ccurve_redshift;
    }
  else if(selection == "ccurveerror")
    {
      return ccurveerror;
    }
  else if(selection == "strongfactorBlk")
    {
      return strongfactorBlk;
    }
  else if(selection == "strongfactorVl")
    {
      return strongfactorVl;
    }
  else if(selection == "datavector")
    {
      return datavector;
    }
  else if(selection == "finaldatavector")
    {
      if(created)
	{
	  return finaldatavector;
	}
    }
  else
    {
      throw invalid_argument("Selection not valid for showdatavector");
    }
}

gsl_vector_int* MPIDataContainer::showintdatavector(string selection)
{
  if(selection == "ccurve")
    {
      return ccurve;
    }
  else
    {
      throw invalid_argument("Selection not valid for showintdatavector");
    }
}
gsl_matrix* MPIDataContainer::showdatamatrix(string selection)
{
  if(selection == "ellip1_covariance")
    {
      return ellip1covariance;
    }
  else if(selection == "ellip2_covariance")
    {
      return ellip2covariance;
    }
  else if(selection == "Fij1")
    {
      return Fij1;
    }
  else if(selection == "Fij2")
    {
      return Fij2;
    }
  else if (selection == "f1_covariance")
    {
      return f1covariance;
    }
  else if(selection == "f2_covariance")
    {
      return f2covariance;
    }
  else if (selection == "g1_covariance")
    {
      return g1covariance;
    }
  else if(selection == "g2_covariance")
    {
      return g2covariance;
    }
  else if(selection == "msysteminfo")
    {
      return msysteminfo;
    }
  else if(selection == "coeff")
    {
      return coeff;
    }
  else if(selection == "finalcoeff")
    {
      if(created)
	{
	  return finalcoeff;
	}
    }
  else
    {
      throw invalid_argument("Selection not valid for showdatamatrix");
    }
}

