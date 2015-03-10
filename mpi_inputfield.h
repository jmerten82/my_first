#ifndef   	INPUTFIELD_H_
# define   	INPUTFIELD_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <cmath>
#include <ctime>
#include "rw_fits.h"
#include "util.h"
#include "soph_math.h"
#include "options.h"
#include "masked_fin_dif.h"
#include "mpi.h"
#include "mpi_comm.h"

using namespace std;

//auxilliary functions

string parametercut2(string);  
// cuts out substring between two kinds of given string

int countsymbol2(string, const char&); 
// counts the given char in a string


void cutints2(string, string, gsl_vector_int*,int);
//gives all integers separated by string in a string to gsl_vector_int 
						
void cutdoubles2(string, string, gsl_vector*,int);
//same for doubles


class DataGrid
{

  /**
     Central class of the field preparation process which contains and
     calculates all relevant quantities which are written to a FITS file later.
  **/


 private:

  double fieldx;
  /* 
     physical x-dimension of the DataGrid.
  */

  double fieldy;
  /* 
     physcial y-dimension of the DataGrid.
  */

  double fieldsize;
  /*
    physical area of the DataGrid, fieldy*fieldy
  */

  int x_dim;
  /*
    pixel x-dimension of the DataGrid
  */

  int y_dim;
  /*
    pixel y-dimension of the DataGrid
  */

  double xcentre;
  /*
    physical x-coordinate of the centre of the bottom left pixel of the DataGrid
  */

  double ycentre;
  /*
    physical y-coordinate of the centre of the bottom left pixel of the DataGrid
  */

  double pixelsize;
  /*
    physical sidelength of one pixel of the DataGrid
  */

  int fieldpixels;
  /*
    Number of non-masked pixels in the DataGrid
  */

  gsl_matrix_int *maskcheck;
  /*
    Identifier for masked pixels in the DataGrid
  */

  bool desshear;
  /*
    For memory deallocation reasons. States if the DataGrid calculates shear.
  */

  bool desflexion;
  /*
    For memory deallocation reasons. States if the DataGrid calculates flexion.
  */
  bool desccurve;
  /*
    For memory deallocation reasons. States if the DataGrid calculates a 
    critical curve estimator.
  */

 
  bool desmsystems;
  /*
    For memory deallocation reasons. States if the DataGrid calculates mutliple
    image systems.
  */

  bool smallindex;
  /*
    Again for deallocation
  */

  bool medindex;
  /*
    see above
  */

  bool root;
  /*
    Am I the root process?
  */


  int numgalshear;
  /*
    Number of galaxies used for shear estimation in the DataGrid.
  */

  double galaxydensityshear;
  /*
    Number density of the galaxies used for shear estimation in the DataGrid.
    numgalshear/fieldsize
  */

  gsl_vector *rec;
  /*
    Contains the RA coordinates of the shear catalogue.
  */

  gsl_vector *dec;
 /*
    Contains the declination coordinates of the shear catalogue.
  */

  gsl_vector *ellip1;
  /*
    Contains the first reduced shear component of the shear catalogue.
  */

  gsl_vector *ellip2;
  /*
    Contains the second reduced shear component of the shear catalogue.
  */

  gsl_vector *shearweight;
  /*
    Contains the weighting function for the reduced shear in the 
    shear catalogue.
  */

  gsl_matrix *meanellip1;
  /* Contains the averaged first reduced shear component for each pixel 
     of the DataGrid.
  */

  gsl_matrix *finalmeanellip1;
  /*
    The full matrix, patched together with all the MPI pieces.
  */

  gsl_matrix *meanellip2;
  /* Contains the averaged second reduced shear component for each pixel 
     of the DataGrid.
  */
  gsl_matrix *finalmeanellip2;
  /*
    The full matrix, patched together with all the MPI pieces.
  */

  gsl_matrix *ellip1sd;
  /* Contains the standard deviation of the first reduced shear component for 
     each pixel of the DataGrid.
  */
  gsl_matrix *finalellip1sd;
  /*
    The full matrix, patched together with all the MPI pieces.
  */

  gsl_matrix *ellip2sd;
  /* Contains the standard deviation of the second reduced shear component for 
     each pixel of the DataGrid.
  */
  gsl_matrix *finalellip2sd;
  /*
    The full matrix, patched together with all the MPI pieces.
  */
  gsl_vector *internalellip1;
  /* Contains the averaged first reduced shear component for each non-masked 
     pixel of the DataGrid.
  */

  gsl_vector *internalellip2;
  /* Contains the averaged second reduced shear component for each non-masked 
     pixel of the DataGrid.
  */
  gsl_vector *internalellip1sd;
  /* Contains the standard deviation of the first reduced shear component for 
     each non-masked pixel of the DataGrid.
  */
  gsl_vector *internalellip2sd;
  /* Contains the standard deviation of the second reduced shear component for 
     each non-masked pixel of the DataGrid.
  */
  gsl_matrix_uchar *galaxyownersellip;
  /*
    Matrix which indicates which galaxy in the original catalogue is used
    for each pixel in the DataGrid. This is needed to perform covariances.
  */
  gsl_matrix_uchar *finalgalaxyownersellip;
  /*
    The gathering matrix which collects from all processes.
  */

  gsl_matrix_int *galaxyshareellip;
  /*
    Matrix which visualises the galaxy overlap in each pixel. This is used
    to calculate the covariances.
  */
  gsl_matrix_int *finalgalaxyshareellip;
  /*
    The full matrix, patched together with all the MPI pieces.
  */

  gsl_matrix *ellip1covariance;
  /*
    Covariance matrix for thr first shear component.
  */

  gsl_matrix *finalellip1covariance;
  /*
    The full matrix, patched together with all the MPI pieces.
  */
  gsl_matrix *ellip2covariance;
  /*
    Covariance matrix for thr second shear component.
  */
  gsl_matrix *finalellip2covariance;
  /*
    The full matrix, patched together with all the MPI pieces.
  */

  int numgalflexion;
  /*
    Number of galaxies used for flexion estimation in the DataGrid.
  */

  double galaxydensityflexion;
  /*
    Number density of the galaxies used for flexion estimation in the DataGrid.
    numgalflexion/fieldsize
  */

  gsl_vector *recflexion;
  /*
    Contains the RA coordinates of the flexion catalogue.
  */

  gsl_vector *decflexion;
  /*
    Contains the Declination coordinates of the flexion catalogue.
  */

  gsl_vector *f1;
  /*
    Contains the first F component of the flexion catalogue.
  */

  gsl_vector *f2;
  /*
    Contains the second F component of the flexion catalogue.
  */

  gsl_vector *g1;
  /*
    Contains the first G component of the flexion catalogue.
  */

  gsl_vector *g2;
  /*
    Contains the second G component of the flexion catalogue.
  */

  gsl_vector *flexionweight;
  /*
    Contains the weighting function for the flexion in the 
    flexion catalogue.
  */

  gsl_matrix *meanf1;
  /* Contains the averaged first F component for each pixel 
     of the DataGrid.
  */
  gsl_matrix *finalmeanf1;

  gsl_matrix *meanf2;
  /* Contains the averaged second F component for each pixel 
     of the DataGrid.
  */
  gsl_matrix *finalmeanf2;


  gsl_matrix *f1sd;
  /* Contains the standard deviation of the first F component for 
     each pixel of the DataGrid.
  */
  gsl_matrix *finalf1sd;


  gsl_matrix *f2sd;
  /* Contains the standard deviation of the second Fcomponent for 
     each pixel of the DataGrid.
  */
  gsl_matrix *finalf2sd;


  gsl_matrix *meang1;
  /* Contains the averaged first G component for each pixel 
     of the DataGrid.
  */
  gsl_matrix *finalmeang1;


  gsl_matrix *meang2;
  /* Contains the averaged second G component for each pixel 
     of the DataGrid.
  */
  gsl_matrix *finalmeang2;


  gsl_matrix *g1sd;
  /* Contains the standard deviation of the first G component for 
     each pixel of the DataGrid.
  */
  gsl_matrix *finalg1sd;

  gsl_matrix *g2sd;
  /* Contains the standard deviation of the second G component for 
     each pixel of the DataGrid.
  */
  gsl_matrix *finalg2sd;

  gsl_vector *internalf1;
  /* Contains the averaged first F component for each non-masked 
     pixel of the DataGrid.
  */
  gsl_vector *internalf2;
  /* Contains the averaged second F component for each non-masked 
     pixel of the DataGrid.
  */
  gsl_vector *internalf1sd;
  /* Contains the standard deviation of the first F component for 
     each non-masked pixel of the DataGrid.
  */
  gsl_vector *internalf2sd;
  /* Contains the standard deviation of the second F component for 
     each non-masked pixel of the DataGrid.
  */
  gsl_vector *internalg1;
  /* Contains the averaged first G component for each non-masked 
     pixel of the DataGrid.
  */
  gsl_vector *internalg2;
  /* Contains the averaged second G component for each non-masked 
     pixel of the DataGrid.
  */
  gsl_vector *internalg1sd;
  /* Contains the standard deviation of the first G component for 
     each non-masked pixel of the DataGrid.
  */
  gsl_vector *internalg2sd;
  /* Contains the standard deviation of the second G component for 
     each non-masked pixel of the DataGrid.
  */
  gsl_matrix_uchar *galaxyownersflexion;
  /*
    Matrix which indicates which galaxy in the original catalogue is used
    for each pixel in the DataGrid. This is needed to perform covariances.
  */
  gsl_matrix_uchar *finalgalaxyownersflexion;
  /*
    Gathering matrix for galaxyowners.
  */

  gsl_matrix_int *galaxyshareflexion;
  /*
    Matrix which visualises the galaxy overlap in each pixel. This is used
    to calculate the covariances.
  */

  gsl_matrix_int *finalgalaxyshareflexion;


  gsl_matrix *f1covariance;
  /*
    Covariance matrix for the first F component.
  */
  gsl_matrix *finalf1covariance;


  gsl_matrix *f2covariance;
  /*
    Covariance matrix for the second F component.
  */
  gsl_matrix *finalf2covariance;


  gsl_matrix *g1covariance;
  /*
    Covariance matrix for the first G component.
  */
  gsl_matrix *finalg1covariance;


  gsl_matrix *g2covariance;
  /*
    Covariance matrix for the second G component.
  */
  gsl_matrix *finalg2covariance;


  int numccurvepts;
  /*
    Number of points in the DataGrid, indicated as part of the critical
    curve estimate.
  */

  gsl_vector *ccurverec;
  /*
    RA coordinates of the crtical curve estimators.
  */

  gsl_vector *ccurvedec;
  /*
    Declination coordinates of the crtical curve estimators.
  */

  gsl_vector *ccurvered;
  /*
    Redshift information on the critical curve estimators.
  */

  gsl_matrix_int *ccurve;
  /*
    Indicator for the actual critical curve estimate position in the
    final DataGrid.
  */

  gsl_matrix *ccurveerror;
  /*
    Error on the critical curve estimation.
  */

  gsl_matrix *ccurveredshift;
  /*
    Redshift information on the critical curve estimator in the actual
    DataGrid.
  */

  int nummsystems;
  /*
    Number of multiple image systems.
  */

  gsl_matrix *msysteminfo;
  /*
    Abstract data matrix containing information about position, gridposition
    and redshift of multiple image systems.
  */

  /*
    Insert for an issue in Charles' project! NOT MEMORY EFFICIENT!
  */

  gsl_matrix *pure_covariance1;
  gsl_matrix *pure_covariance2;


 public:

  DataGrid(FieldOptions&,int,int my_rank,int p);
  /*
    Constructor, which needs an options class and an iterationindex.
    Reads the catalogues and allocates memory, the rank makes sure
    that only the leading process reads data.
  */

  ~DataGrid();
  /*
    Standard destructor, free (a lot of) memory.
  */

  void analyse(FieldOptions&,int my_rank, int p);
  /*
    Performs averaging pixel associations and covariances.
    The time consuming part.  
  */

  void write(FieldOptions&,int,int my_rank);
  /*
    Write the results to FITS resp. ASCII, depending on the options
    class and the iterationindex which contain the right filenames.
  */

};


#endif 	    /* !INPUTFIELD_H_ */
