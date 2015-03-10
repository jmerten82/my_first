#include "mpi_coreroutines.h"

using namespace std;

void MPIWeakFactor(ReconstructionOptions &options,GalaxyCluster &cluster, MPIDataContainer &container)
{
  if(options.showbool("shear"))
    {
      double shearz = cosmicweight(options.showdouble("clusterz"),gsl_vector_get(container.showdatavector("shear_redshift"),0));
      for(int i = 0; i < container.showdatamatrix("Fij1")->size1; i++)
	{
	  for(int j = 0; j < container.showdatamatrix("Fij1")->size2; j++)
	    {
	      gsl_matrix_set(container.showdatamatrix("Fij1"),i,j,shearz*shearz*gsl_matrix_get(container.showdatamatrix("ellip1_covariance"),i,j)*1.0/((1.0-shearz*gsl_vector_get(cluster.data("convergence"),i))*(1.0-shearz*gsl_vector_get(cluster.data("convergence"),j))));
	      gsl_matrix_set(container.showdatamatrix("Fij2"),i,j,shearz*shearz*gsl_matrix_get(container.showdatamatrix("ellip2_covariance"),i,j)*1.0/((1.0-shearz*gsl_vector_get(cluster.data("convergence"),i))*(1.0-shearz*gsl_vector_get(cluster.data("convergence"),j))));
	    }
	}
    }
  if(options.showbool("flexion"))
    {
      double flexionz = cosmicweight(options.showdouble("clusterz"),gsl_vector_get(container.showdatavector("flexion_redshift"),0));

      for(int i = 0; i < container.showdatamatrix("f1_covariance")->size1; i++)
	{
	  for(int j = 0; j < container.showdatamatrix("f1_covariance")->size2; j++)
	    {
	      gsl_matrix_set(container.showdatamatrix("f1_covariance"),i,j,gsl_matrix_get(container.showdatamatrix("f1_covariance"),i,j)*flexionz*flexionz);
	      gsl_matrix_set(container.showdatamatrix("f2_covariance"),i,j,gsl_matrix_get(container.showdatamatrix("f2_covariance"),i,j)*flexionz*flexionz);
	      gsl_matrix_set(container.showdatamatrix("g1_covariance"),i,j,gsl_matrix_get(container.showdatamatrix("g1_covariance"),i,j)*flexionz*flexionz);
	      gsl_matrix_set(container.showdatamatrix("g2_covariance"),i,j,gsl_matrix_get(container.showdatamatrix("g2_covariance"),i,j)*flexionz*flexionz);
	    }
	}
    }
}

void MPIStrongFactor(ReconstructionOptions &options,GalaxyCluster &cluster, MPIDataContainer &container)
{
  if(options.showbool("ccurve"))
    {
      double error = gsl_vector_get(container.showdatavector("ccurveerror"),0);

      for(int i = 0; i < container.showdatavector("strongfactorBlk")->size; i++)
	{
	  if(gsl_vector_int_get(container.showintdatavector("ccurve"),i) != 0)
	    {
	      gsl_vector_set(container.showdatavector("strongfactorBlk"),i,4.0*gsl_vector_get(cluster.data("jacdet"),i)*pow(cosmicweight(options.showdouble("clusterz"),gsl_vector_get(container.showdatavector("ccurve_redshift"),i)),2)*1/pow(error,2));
	      gsl_vector_set(container.showdatavector("strongfactorVl"),i,4.0*gsl_vector_get(cluster.data("jacdet"),i)*cosmicweight(options.showdouble("clusterz"),gsl_vector_get(container.showdatavector("ccurve_redshift"),i))*1/pow(error,2));
	    }
	}
    }
}


void MPI_Blk(ReconstructionOptions &options,GalaxyCluster &referencecluster,MPIDataContainer &container,int p, int my_rank)
{
  FinDifGrid recgrid(referencecluster.grid(),referencecluster.showdouble("x_frct"),true);

  MPI_Array_Manager array1(my_rank,p,recgrid.showint("fieldpixels"),0);

  if(my_rank == 0)
    {
      cout <<"Building up coefficient matrix" <<endl;
      cout <<"Progress (of process 0):   |         |" <<endl;
      CURSORUP(1);
      CURSORRIGHT(28);
    }
  double value = 0.0;
  double progress = 0.1;

  for(int l = 0; l < array1.showsize(); l++)
    {
      int ll = l + array1.showstartindex();
      if(l < recgrid.showint("fieldpixels"))
	{
	  if(my_rank == 0)
	    {
	      if( (double) l/array1.showsize() > progress)
		{
		  cout <<"#" <<flush;
		  progress += 0.1;
		}
	    }
	  for(int k = 0; k < recgrid.showint("fieldpixels"); k++)
	    {  
	      value = 0.0;
	      if(options.showbool("shear"))
		{
		    if(options.showbool("weaklimit"))
			{
			    value += recgrid.ShearBlk(ll,k,container.showdatavector("ellip1"),container.showdatavector("ellip2"),container.showdatamatrix("ellip1_covariance"),container.showdatamatrix("ellip2_covariance"),"shear");
			}
		    else
			{
			    value += recgrid.ShearBlk(ll,k,container.showdatavector("ellip1"),container.showdatavector("ellip2"),container.showdatamatrix("Fij1"),container.showdatamatrix("Fij2"),"reduced shear");
			}
		}
	      if(options.showbool("flexion"))
		{
		  value += recgrid.FlexionBlk(ll,k,container.showdatamatrix("f1_covariance"),container.showdatamatrix("f2_covariance"),container.showdatamatrix("g1_covariance"),container.showdatamatrix("g2_covariance"));
		}
	      if(options.showbool("ccurve"))
		{
		  value += recgrid.CcurveBlk(ll,k,container.showintdatavector("ccurve"),container.showdatavector("strongfactorBlk"));
		}
	      if(options.showbool("msystems"))
		{
		  value += recgrid.MsystemBlk(ll,k,container.showdatamatrix("msysteminfo"));
		}
	      if(recgrid.showint("x_dim") == options.showint("startdim"))
		{
		  value += recgrid.RegularisationBlk(ll,k,options.showdouble("regshear"),options.showdouble("regflexion"),"convergence");		  
		}
	      else
		{
		  value += recgrid.RegularisationBlk(ll,k,options.showdouble("regshear"),options.showdouble("regflexion"),options.showstring("regtype"));
		}

	      gsl_matrix_set(container.showdatamatrix("coeff"),l,k,value);
	    }
	}
    }
	
  if(my_rank == 0)
    {
      cout <<endl;
    }
}


void MPI_Vl(ReconstructionOptions &options,GalaxyCluster &referencecluster,MPIDataContainer &container,int p, int my_rank)
{
  FinDifGrid recgrid(referencecluster.grid(),referencecluster.showdouble("x_frct"),true);
  MPI_Array_Manager array1(my_rank,p,recgrid.showint("fieldpixels"),0);  

  if(my_rank == 0)
    {  
      cout <<"Building up data vector" <<endl;
      cout <<"Progress (of slowest process):   |         |" <<endl;
      CURSORUP(1);
      CURSORRIGHT(34);
    }
  double value = 0.0;
  double progress = 0.1;
  

  
  for(int l = 0; l < array1.showsize(); l++)
    {
      int ll = l + array1.showstartindex();
      if(l < recgrid.showint("fieldpixels"))
	{
	  if(my_rank == 0)
	    {
	      if( (double) l/array1.showsize() > progress)
		{
		  cout <<"#" <<flush;
		  progress += 0.1;
		}
	    }      
	  value = 0.0;
	  if(options.showbool("shear"))
	    {
		if(options.showbool("weaklimit"))
		    {
			value += recgrid.ShearVl(ll,container.showdatavector("ellip1"),container.showdatavector("ellip2"),container.showdatamatrix("ellip1_covariance"),container.showdatamatrix("ellip2_covariance"),container.showdatavector("shear_redshift"),"shear");
		    }
		else
		    {
			value += recgrid.ShearVl(ll,container.showdatavector("ellip1"),container.showdatavector("ellip2"),container.showdatamatrix("Fij1"),container.showdatamatrix("Fij2"),container.showdatavector("shear_redshift"),"reduced shear");
		    }
	    }
	  if(options.showbool("flexion"))
	    {
	      value += recgrid.FlexionVl(ll,container.showdatavector("f1"),container.showdatavector("f2"),container.showdatavector("f1"),container.showdatavector("g2"),container.showdatamatrix("f1_covariance"),container.showdatamatrix("f2_covariance"),container.showdatamatrix("g1_covariance"),container.showdatamatrix("g2_covariance"),container.showdatavector("flexion_redshift"));
	    }
	  if(options.showbool("ccurve"))
	    {
	      value += recgrid.CcurveVl(ll,container.showintdatavector("ccurve"),container.showdatavector("strongfactorVl"));
	    }
	  if(options.showbool("msystems"))
	    {
	      value += recgrid.MsystemVl(ll,container.showdatamatrix("msysteminfo"));
	    }
	  if(recgrid.showint("x_dim") == options.showint("startdim"))
	    {
	      value += recgrid.RegularisationVl(ll,referencecluster.data("convergence"),referencecluster.data("shear1"),referencecluster.data("shear2"),referencecluster.data("f1"),referencecluster.data("f2"),referencecluster.data("g1"),referencecluster.data("g2"),options.showdouble("regshear"),options.showdouble("regflexion"),"convergence");
	    }
	  else
	    {
	      value += recgrid.RegularisationVl(ll,referencecluster.data("convergence"),referencecluster.data("shear1"),referencecluster.data("shear2"),referencecluster.data("f1"),referencecluster.data("f2"),referencecluster.data("g1"),referencecluster.data("g2"),options.showdouble("regshear"),options.showdouble("regflexion"),options.showstring("regtype"));
	    }

	  gsl_vector_set(container.showdatavector("datavector"),l,value);
	}
    }

  if(my_rank == 0)
    {
      cout <<endl;
    }
}








  











