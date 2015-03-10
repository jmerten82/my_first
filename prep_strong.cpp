#include "prep_strong.h"

using namespace std;


void masked_interpolcut(HighresOptions &options1)
{

  int ox_dim = read_intheader(options1.show_filename("lowres_in"),"NAXIS1");
  int oy_dim = read_intheader(options1.show_filename("lowres_in"),"NAXIS2");
  int ix_dim = options1.showint("finalx");
  int iy_dim = options1.showint("finaly");
  int cutx = options1.showint("x_max")-options1.showint("x_min");
  int cuty = options1.showint("y_max")-options1.showint("y_min");
  int zoom_dim = cutx*cuty;


  gsl_matrix *inpot = gsl_matrix_calloc(oy_dim,ox_dim);
  gsl_matrix_int *inmask = gsl_matrix_int_calloc(oy_dim,ox_dim);
  gsl_matrix *bigpot = gsl_matrix_calloc(iy_dim,ix_dim);
  gsl_matrix *outpot = gsl_matrix_calloc(cuty,cutx);
  gsl_matrix_int *bigccurve;
  gsl_matrix *bigredshift;
  gsl_matrix_int *outccurve;
  gsl_matrix *outredshift;
  gsl_vector_int *numbers;
  gsl_matrix *msystems;
  gsl_vector_int *cutnumbers;
  gsl_matrix_int *cutnumbers2;


  read_pimg(options1.show_filename("lowres_in"),inpot);
  read_imgeint(options1.show_filename("lowres_in"),"field_mask",inmask);

  cout <<"Interpolating potential..." <<flush;
  cspline_smooth_morecleverinterpol(inpot,inmask,bigpot);
  cout <<"Done." <<endl;

  if(options1.showbool("ccurve"))
    {
      bigccurve = gsl_matrix_int_calloc(iy_dim,ix_dim);
      bigredshift = gsl_matrix_calloc(iy_dim,ix_dim);
      outccurve = gsl_matrix_int_calloc(cuty,cutx);
      outredshift = gsl_matrix_calloc(cuty,cutx);
      read_pimgint(options1.show_filename("ccurve_in"),bigccurve);
      read_imge(options1.show_filename("ccurve_in"),"redshift_info",bigredshift);
    }

  if(options1.showbool("msystems"))
    {
      int i_dim = options1.showint("finalx")*options1.showint("finaly");
      
      numbers = gsl_vector_int_calloc(i_dim);
      cutnumbers = gsl_vector_int_calloc(zoom_dim);
      for(int i = 0; i < i_dim; i++)
	{
	  gsl_vector_int_set(numbers,i,i);
	}
      
      msystems = gsl_matrix_calloc(18,options1.showint("nummsystems"));
      ReadMsystemInfo(options1.show_filename("msystem_in"),msystems);
    }

  cout <<"Cutting files..." <<flush;

  cut(bigpot,ix_dim,iy_dim,options1.showint("x_min"),options1.showint("x_max"),options1.showint("y_min"),options1.showint("y_max"),outpot);

  if(options1.showbool("ccurve"))
    {
      cut(bigccurve,ix_dim,iy_dim,options1.showint("x_min"),options1.showint("x_max"),options1.showint("y_min"),options1.showint("y_max"),outccurve);
      cut(bigredshift,ix_dim,iy_dim,options1.showint("x_min"),options1.showint("x_max"),options1.showint("y_min"),options1.showint("y_max"),outredshift);
    }
  if(options1.showbool("msystems"))
    {
      cut(numbers,ix_dim,iy_dim,options1.showint("x_min"),options1.showint("x_max"),options1.showint("y_min"),options1.showint("y_max"),cutnumbers);

      cutnumbers2 = gsl_matrix_int_calloc(cuty,cutx);
      vecmat(cutnumbers,cutnumbers2);
      
      int counter = 0;
      for(int ms = 0; ms < options1.showint("nummsystems"); ms++)
	{
	  for(int s = 0; s < (int) gsl_matrix_get(msystems,0,ms); s++)
	    {
	      for(int i = 0; i < cuty; i++)
		{
		  for(int j = 0;  j < cutx; j++)
		    {
		      if(gsl_matrix_int_get(cutnumbers2,i,j) == (int) gsl_matrix_get(msystems,12+s,ms))
			{
			  gsl_matrix_set(msystems,12+s,ms,counter);
			}
		      counter++;
		    }
		}
	      counter = 0;
	    }
	}
    }
  cout <<"Done." <<endl;

  cout <<"Writing files..." <<flush;

  write_pimg(options1.show_filename("core_out"),bigpot);
  write_imge(options1.show_filename("core_out"),"cut_potential",outpot);

  if(options1.showbool("ccurve"))
    {
      write_imgeint(options1.show_filename("core_out"),"cut_ccurve",outccurve);
      write_imge(options1.show_filename("core_out"),"cut_redshift",outredshift);
    }
  if(options1.showbool("msystems"))
    {
	WriteMsystemInfo(msystems,options1.show_filename("msystem_reset"));
    }
  cout <<"Done." <<endl;
}
