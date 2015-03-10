#include <iostream>
#include <tclap/CmdLine.h>
#include <gsl/gsl_matrix.h>
#include <rw_fits.h>
#include <util.h>

using namespace std;

int main(int argc, char* argv[])
{

  TCLAP::CmdLine cmd("Cosmic weight calculation", ' ', "0.1");
  TCLAP::ValueArg<std::string> inputArg("i","input","Input redshift map",true,"","string",cmd);
  TCLAP::ValueArg<std::string> outputArg("o","output","Output cosmological weight map",true,"","string",cmd);
  TCLAP::ValueArg<double> zArg("z","z_cluster","Lens redshift",true,0.0,"double",cmd);
  TCLAP::ValueArg<double> sArg("s","s_sources","Weak lensing source redshift",true,0.0,"double",cmd);

  cmd.parse( argc, argv );

  string infile = inputArg.getValue();
  string outfile = outputArg.getValue();
  double zlens = zArg.getValue();
  double zsource = sArg.getValue();

  double x_dim = read_intheader(infile,"NAXIS1");
  double y_dim = read_intheader(infile,"NAXIS2");

  gsl_matrix *zmap = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *cmap = gsl_matrix_calloc(y_dim,x_dim);

  double csource = cosmicweight(zlens,zsource);

  read_imge(infile,"redshift_info",zmap);

  gsl_matrix_set_all(cmap,csource);

  for(int i = 0; i < y_dim; i++)
    {
      for(int j = 0; j < x_dim; j++)
	{
	  if(gsl_matrix_get(zmap,i,j) != 0.0)
	    {
	      gsl_matrix_set(cmap,i,j,cosmicweight(zlens,gsl_matrix_get(zmap,i,j)));
	    }
	  else
	    {
	      gsl_matrix_set(zmap,i,j,zsource);
	    }
	}
    }

  write_pimg(outfile,zmap);
  write_imge(outfile,"cosmic_weight",cmap);

  return 0;

}
