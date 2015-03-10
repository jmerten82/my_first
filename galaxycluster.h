#ifndef   	GALAXYCLUSTER_H
# define   	GALAXYCLUSTER_H

#include <stdexcept>
#include <gsl/gsl_matrix.h>
#include <string>
#include <ctime>
#include "masked_fin_dif.h"
#include "util.h"
#include "rw_fits.h"
#include "options.h"


using namespace std;

class GalaxyCluster
{

  /**
     Represents in general a lens with certain lensing properties,
     these are calculated from a lensing potential via a finite differencing
     approach. The whole construct lives on a well-defined grid, which can also
     be masked.
   **/

 private:

  int x_dim;
  /* 
     x-dimension of the grid
  */

  int y_dim; 
  /*
    y-dimension of the grid
  */

  double x_frct;
  /* if the lens grid is just a cut-out of a bigger grid, this is the 
     fraction of it's x-dimension to the one of the big field.
   */
 
  int fieldpixels;
  /* number of pixels in the lens grid which are not masked and can be used
     in a reconstruction.
  */
 
  gsl_matrix_int *fieldmask;
  /* vector defining the masking of the field. A value of 1 means that a
a pixel is masked, anything else that it is not masked.
  */

  gsl_vector *pot; 
  /*
    lensing potential
  */

  gsl_vector *a1; 
  /*
    first component of the deflection angle
  */

  gsl_vector *a2;
  /*
    second component of the deflection angle
  */ 
  gsl_vector *shear1;
  /*
    first component of the shear
  */
 
  gsl_vector *shear2; 
  /*
    second of the shear
  */

  gsl_vector *convergence;
  /*
    convergence of the lens, not necessarily normalised
  */

  gsl_vector *f1; 
  /*
    F1 component of flexion
  */

  gsl_vector *f2;
  /*
    F2 component of flexion
  */

  gsl_vector *g1; 
  /*
    G1 component of flexion
  */

  gsl_vector *g2;
  /*
    G2 component of flexion 
  */

  gsl_vector *jacdet; 
  /*
    determinant of the lensing Jacobian
  */

 public:

  GalaxyCluster(gsl_matrix_int *fieldmask,double x_frct);
  /* Constructor, needs a fieldmask in format: value 1 means masked else means
     not masked and x_frct as the fraction of the lens field to the whole
     reconstruction field.
  */

  ~GalaxyCluster();
  /* Standard destructor which deallocates the gsl_vectors
   */

  int showint(string selection);
  /* return the integer cluster properties depending on the selection parameter.
     Selections are "x_dim", "y_dim" and "dim" (x_dim*y_dim)
  */

  double showdouble(string selection);
  /*
    returns the double cluster properties depending on the selection.
    Selection right now is just "x_frct"
  */

  gsl_vector* data(string selection);
  /*
    returns a pointer to the internal data structures of the cluster.
    Selections are: "pot","shear1/2", "convergence", "f1/2", "g1/2" and "jacdet"
  */
 
  gsl_matrix_int* grid();
  /*
    returns a point to the fieldmask of the cluster.
  */

  void buildfrompot();
  /* 
     calculates all the other quantities from a lensing potential
     which was already read by the cluster via readgridvector or readmap
  */

  void buildwithoutpot();
  /* 
     computes the jacdet out of conv and shear
  */

  void masssheetnormalise(string selection); //selection is zero
  /*
    normalises the shear and convergence of the cluster, only selection by
    now is "zero" which means that the pixel with the lowest convergence value
    is set to zero and the rest is normalised accordingly.
  */

  void readgridvector(gsl_vector* input, string selection);
  /* 
     reads an external gsl_vector input into the selection which are:
     "pot","shear1/2", "convergence", "f1/2", "g1/2" and "jacdet".
     Length of these vectors must be fieldpixels.
  */

  void readmap(gsl_matrix* input,string selection);
  /* 
     reads a gsl_matrix into selection which are:
     "pot","shear1/2", "convergence", "f1/2", "g1/2" and "jacdet".
     Row size of this matrix must be y_dim and col size must be x_dim
  */
     
  void writegridvector(gsl_vector* output,string selection);
  /*
    writes the selection into an external gsl_vector output 
    of length fieldpixels. Selections are "pot","shear1/2", 
    "convergence", "f1/2", "g1/2" and "jacdet".
  */

  void writemap(gsl_matrix *output,string selection);
  /*
    writes the selection into an external gsl_matrix output of
    row size y_dim and col size x_dim. Selections are "pot","shear1/2", 
    "convergence", "f1/2", "g1/2" and "jacdet".
  */

  void gridvectortomap(gsl_vector *input,gsl_matrix *output);
  /*
    Converts an external gsl_vector input of length fieldpixels into an
    external gsl_matrix output of row size y_dim and col size x_dim.
  */
 
  void maptogridvector(gsl_matrix* input,gsl_vector *output);
  /*
    Converts an external gsl_matrix of row size y_dim and col size y_dim
    into an external gsl_vector of length fieldpixels.
  */
 
  void gridvectortomap(gsl_vector_int*,gsl_matrix_int*);
  /*
    Converts an external gsl_vector input of length fieldpixels into an
    external gsl_matrix output of row size y_dim and col size x_dim.
  */

  void maptogridvector(gsl_matrix_int*,gsl_vector_int*);
  /*
    Converts an external gsl_matrix of row size y_dim and col size y_dim
    into an external gsl_vector of length fieldpixels.
  */

  void writetofits(string filename);
  /* 
     Writes all the data structures of the cluster into an multi extension
     FITS file with name filename. Masked areas are set to a fixed value
     in the maps, usually 0.0
  */

  void writetofits(ReconstructionOptions &options, int iterationindex);
  /* 
     Writes all the data structures of the cluster into an multi extension
     FITS file with name filename. Masked areas are set to a fixed value
     in the maps, usually 0.0.
     Additionally, this routine writes information about the parameters of
     the reconstruction from which the cluster was derived into the FITS header
     of the pHDU.
  */

};

#endif 	    /* !GALAXYCLUSTER_H */
