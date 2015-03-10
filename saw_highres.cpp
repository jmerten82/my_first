#include <iostream>
#include <cmath>
#include <string>
#include <ctime>
#include <gsl/gsl_matrix.h>
#include "rw_fits.h"
#include "galaxycluster.h"
#include "reconstruction_kernel.h"
#include "prep_strong.h"
#include "options.h"
#include "util.h"
#include "soph_math.h"

using namespace std;

int main(int argc, char *argv[])
{

  cout <<"----------SAW 1.4 High resolution reconstruction routine----------" <<endl;
  cout <<endl;

  //reading in the options file

  HighresOptions options1(argv[1]);

  //defining necessary stuff for the interpolation of the initial field

  int cutx = options1.showint("x_max")-options1.showint("x_min");
  int cuty = options1.showint("y_max")-options1.showint("y_min");
  int dim = cutx*cuty;

  int ox_dim = read_intheader(options1.show_filename("lowres_in"),"NAXIS1");
  int oy_dim = read_intheader(options1.show_filename("lowres_in"),"NAXIS2");

  double frct = (double) cutx/options1.showint("finalx");
  gsl_vector_int *zoomccurve = gsl_vector_int_calloc(dim);
  gsl_vector *zoomredshift = gsl_vector_calloc(dim);
  gsl_matrix *msystemsinfo = gsl_matrix_calloc(18,options1.showint("nummsystems"));
  gsl_vector *pot = gsl_vector_calloc(dim);
  gsl_vector *strongfactorBlk = gsl_vector_calloc(dim);
  gsl_vector *strongfactorVl = gsl_vector_calloc(dim);
  gsl_matrix_int *fakemask = gsl_matrix_int_calloc(cuty,cutx);
  gsl_matrix *Blk = gsl_matrix_calloc(dim,dim);
  gsl_vector *Vl = gsl_vector_calloc(dim);


  //Performing the interpolation, cutting and masking of the initial field
  cout <<"Done" <<endl;
  masked_interpolcut(options1);

  //Reading data for the final high-resolution reconstruction

  read_imge(options1.show_filename("core_out"),"cut_potential",pot);
  if(options1.showbool("ccurve"))
    {
      read_imgeint(options1.show_filename("core_out"),"cut_ccurve",zoomccurve);
      read_imge(options1.show_filename("core_out"),"cut_redshift",zoomredshift);
    }
  if(options1.showbool("msystems"))
    {
      ReadMsystemInfo(options1.show_filename("msystem_reset"),msystemsinfo);
    }

  //Building a reference from the interpolated result
  cout <<frct <<endl;
  GalaxyCluster reference(fakemask,frct);
  reference.readgridvector(pot,"pot");
  reference.buildfrompot();
  //reference.writetofits("test.fits");

  //Performing the reconstruction
  StrongFactor(reference,zoomredshift,options1.showdouble("ssigma"),strongfactorBlk,strongfactorVl);
  HighresBlk(options1,frct,zoomccurve,msystemsinfo,strongfactorBlk,Blk);
  HighresVl(options1,reference,frct,zoomccurve,msystemsinfo,strongfactorVl,Vl);
  cout <<endl;
  cout <<"Solving linear system..." <<flush;
  solve_gsl(Blk,Vl);
  cout <<"Done." <<endl;

  reference.readgridvector(Vl,"pot");
  reference.buildfrompot();

  //Merging the results into a single map with different effective pixel
  //resolutions

  cout <<"Merging results..." <<flush;
  gsl_matrix *convmerge = gsl_matrix_calloc(options1.showint("finaly"),options1.showint("finalx"));
  gsl_matrix *lowconv = gsl_matrix_calloc(oy_dim,ox_dim);
  gsl_matrix *highconv = gsl_matrix_calloc(cuty,cutx);
  gsl_matrix *shear1merge = gsl_matrix_calloc(options1.showint("finaly"),options1.showint("finalx"));
  gsl_matrix *lowshear1 = gsl_matrix_calloc(oy_dim,ox_dim);
  gsl_matrix *highshear1 = gsl_matrix_calloc(cuty,cutx);
  gsl_matrix *shear2merge = gsl_matrix_calloc(options1.showint("finaly"),options1.showint("finalx"));
  gsl_matrix *lowshear2 = gsl_matrix_calloc(oy_dim,ox_dim);
  gsl_matrix *highshear2 = gsl_matrix_calloc(cuty,cutx);
  gsl_matrix *jacdetmerge = gsl_matrix_calloc(options1.showint("finaly"),options1.showint("finalx"));
  gsl_matrix *lowjacdet = gsl_matrix_calloc(oy_dim,ox_dim);
  gsl_matrix *highjacdet = gsl_matrix_calloc(cuty,cutx);
  read_imge(options1.show_filename("lowres_in"),"convergence",lowconv);
  reference.writemap(highconv,"convergence");
  //read_imge(options1.showstring("rectoutput"),"convergence",highconv);
  merge_maps(lowconv,highconv,options1.showint("x_min"),options1.showint("x_max"),options1.showint("y_min"),options1.showint("y_max"),convmerge);
  read_imge(options1.show_filename("lowres_in"),"shear1",lowshear1);
  reference.writemap(highshear1,"shear1");
  //read_imge(options1.showstring("rectoutput"),"shear1",highshear1);
  merge_maps(lowshear1,highshear1,options1.showint("x_min"),options1.showint("x_max"),options1.showint("y_min"),options1.showint("y_max"),shear1merge);
  read_imge(options1.show_filename("lowres_in"),"shear2",lowshear2);
  reference.writemap(highshear2,"shear2");
  //read_imge(options1.showstring("rectoutput"),"shear2",highshear2);
  merge_maps(lowshear2,highshear2,options1.showint("x_min"),options1.showint("x_max"),options1.showint("y_min"),options1.showint("y_max"),shear2merge);
  read_imge(options1.show_filename("lowres_in"),"jacdet",lowjacdet);
  reference.writemap(highjacdet,"jacdet");
  //read_imge(options1.showstring("rectoutput"),"jacdet",highjacdet);
  merge_maps(lowjacdet,highjacdet,options1.showint("x_min"),options1.showint("x_max"),options1.showint("y_min"),options1.showint("y_max"),jacdetmerge);

  write_pimg(options1.show_filename("merge_out"),convmerge);
  write_imge(options1.show_filename("merge_out"),"shear1",shear1merge);
  write_imge(options1.show_filename("merge_out"),"shear2",shear2merge);
  write_imge(options1.show_filename("merge_out"),"jacdet",jacdetmerge);
  string method = " ";
  if(options1.showbool("ccurve"))
    {
      method += "c ";
    }
  if(options1.showbool("msystems"))
    {
      method += "m ";
    }

  write_header(options1.show_filename("merge_out"),"METHOD",method," ");

  write_header(options1.show_filename("merge_out"),"LOWRES",options1.showstring("lowresresult")," ");
  write_header(options1.show_filename("merge_out"),"CCURVE",options1.showstring("highresccurve")," ");
  write_header(options1.show_filename("merge_out"),"MSYSTEMS",options1.showstring("highresmsystemsfile")," ");
  write_header(options1.show_filename("merge_out"),"CUTX1",options1.showint("x_min")," ");
  write_header(options1.show_filename("merge_out"),"CUTX2",options1.showint("x_max")," ");
  write_header(options1.show_filename("merge_out"),"CUTY1",options1.showint("y_min")," ");
  write_header(options1.show_filename("merge_out"),"CUTY2",options1.showint("y_max")," ");
  write_header(options1.show_filename("merge_out"),"REG",options1.showdouble("reg")," ");
  write_header(options1.show_filename("merge_out"),"REGTYPE",options1.showstring("regscheme")," ");
  write_header(options1.show_filename("merge_out"),"SIGMA",options1.showdouble("ssigma")," ");
	 
  cout <<"Done." <<endl;
  return 0;
}
