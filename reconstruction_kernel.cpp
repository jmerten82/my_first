#include "reconstruction_kernel.h"

void StrongFactor(GalaxyCluster &cluster,gsl_vector*redshift,double sigma, gsl_vector *strongfactorBlk, gsl_vector *strongfactorVl)
{

  for(int i = 0; i < strongfactorBlk->size; i++)
    {
      gsl_vector_set(strongfactorBlk,i,4.0*gsl_vector_get(cluster.data("jacdet"),i)*gsl_vector_get(redshift,i)*gsl_vector_get(redshift,i)/(sigma*sigma));
      gsl_vector_set(strongfactorVl,i,4.0*gsl_vector_get(cluster.data("jacdet"),i)*gsl_vector_get(redshift,i)*gsl_vector_get(redshift,i)/(sigma*sigma));
    }
}


void HighresBlk(HighresOptions &options1,double frct,gsl_vector_int *ccurve, gsl_matrix *MSystemInfo,gsl_vector *strongfactor,gsl_matrix *Blk)
{
  int cutx = options1.showint("x_max")-options1.showint("x_min");
  int cuty = options1.showint("y_max")-options1.showint("y_min");
  gsl_matrix_int *fakemask = gsl_matrix_int_calloc(cuty,cutx);
  FinDifGrid recgrid (fakemask,frct,true);

  cout <<endl;
  cout <<"Building up coefficient matrix" <<endl;
  cout <<"Progress:   |         |" <<endl;
  CURSORUP(1);
  CURSORRIGHT(13);
  double progress = 0.1;
  double value;

  for(int l = 0; l < Blk->size1; l++)
    {
      if( (double) l/Blk->size1 > progress)
	{
	  cout <<"#" <<flush;
	  progress += 0.1;
	}
      for(int k = 0; k < Blk->size2; k++)
	{
	  value = 0.0;

	  if(options1.showbool("ccurve"))
	    {
	      value += recgrid.CcurveBlk(l,k,ccurve,strongfactor);
	    }
	  if(options1.showbool("msystems"))
	    {
	      value += recgrid.MsystemBlk(l,k,MSystemInfo);
	    }
	  value += recgrid.RegularisationBlk(l,k,options1.showdouble("reg"),options1.showdouble("reg"),options1.showstring("regscheme"));

	  gsl_matrix_set(Blk,l,k,value);
	}
    }
}
void HighresVl(HighresOptions &options1,GalaxyCluster &reference,double frct,gsl_vector_int *ccurve, gsl_matrix *MSystemInfo,gsl_vector *strongfactor,gsl_vector *Vl)
{
  int cutx = options1.showint("x_max")-options1.showint("x_min");
  int cuty = options1.showint("y_max")-options1.showint("y_min");

  gsl_matrix_int *fakemask = gsl_matrix_int_calloc(cuty,cutx);
  FinDifGrid recgrid (fakemask,frct,true);

  cout <<endl;
  cout <<"Building up data vector" <<endl;
  cout <<"Progress:   |         |" <<endl;
  CURSORUP(1);
  CURSORRIGHT(13);
  double progress = 0.1;
  double value;

  for(int l = 0; l < Vl->size; l++)
    {
      if( (double) l/Vl->size > progress)
	{
	  cout <<"#" <<flush;
	  progress += 0.1;
	}
      
      value = 0.0;
      
      if(options1.showbool("ccurve"))
	{
	  value += recgrid.CcurveVl(l,ccurve,strongfactor);
	}
      if(options1.showbool("msystems"))
	{
	  value += recgrid.MsystemVl(l,MSystemInfo);
	}
      value += recgrid.RegularisationVl(l,reference.data("convergence"),reference.data("shear1"),reference.data("shear2"),reference.data("f1"),reference.data("f2"),reference.data("g1"),reference.data("g2"),options1.showdouble("reg"),options1.showdouble("reg"),options1.showstring("regscheme"));
      gsl_vector_set(Vl,l,value);
    }
}

 







