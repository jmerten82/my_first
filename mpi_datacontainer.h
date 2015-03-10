#ifndef   	MPI_DATACONTAINER_H_
# define   	MPI_DATACONTAINER_H_

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <gsl/gsl_matrix.h>
#include "options.h"
#include "mpi.h"
#include "rw_fits.h"
#include "mpi_comm.h"
#include "masked_fin_dif.h"
#include "util.h"

using namespace std;


class MPIDataContainer
{
  /**
     Class which contains all the necessary input data for a reconstruction.
     Since the code uses an MPI driver the class also contains routines
     to send around this data between the processes and reads them only
     on one process.
  **/


 private:

  int dim;
  /*
    The general dimension of the container.
  */

  bool shear;
  /*
    Flag which states if class contains shear data.
  */

  bool flexion;
  /*
    Flag which states if class contains flexion data.
  */
  bool ccurveflag;
  /*
    Flag which states if class contains ccurve data.
  */

  bool msystems;
  /*
    Flag which states if class contains multiple image system data.
  */

  bool created;
  /*
    Takes the process part in matrix build-up?
  */

  bool root;
  /*
    Am I the leading process;
  */

  int my_rank;
  /*
    Global rank of the process.
  */
  int p;
  /*
    Total number of processes.
  */

  gsl_matrix *coeff;
  /*
    The coefficient matrix used in the reconstruction.
  */

  gsl_matrix *finalcoeff;
  /*
    The full coefficient matrix summing all contributions.
  */

  gsl_vector *datavector;
  /*
    The data vector used in the reconstruction.
  */
  gsl_vector *finaldatavector;
  /*
    see above
  */

  gsl_vector *ellip1;
  /*
    Ellipticity component 1 .
  */

  gsl_vector *ellip2;
  /*
    Ellipticity component 1 .
  */

  gsl_matrix *ellip1covariance;
  /* 
     Covariance matrix for ellipticity component 1.
  */

  gsl_matrix *ellip2covariance;
  /* 
     Covariance matrix for ellipticity component 2.
  */

  gsl_matrix *Fij1;
  /*
    Prefactor for the reduced shear 1 which changes in each inner-level
    iteration.
  */

  gsl_matrix *Fij2;
  /*
    Prefactor for the reduced shear 2 which changes in each inner-level
    iteration.
  */

  gsl_vector *shear_redshift;
  /* 
     Shear redshift vector.
  */

  gsl_vector *f1;
  /* 
     Flexion component F1.
  */

  gsl_vector *f2;
  /* 
     Flexion component F2.
  */

  gsl_vector *g1;
  /* 
     Flexion component G1.
  */

  gsl_vector *g2;
  /* 
     Flexion component G2.
  */

  gsl_matrix *f1covariance;
  /* 
     Covariance matrix for flexion component F1.
  */
  gsl_matrix *f2covariance;
  /* 
     Covariance matrix for flexion component F2.
  */

  gsl_matrix *g1covariance;
  /* 
     Covariance matrix for flexion component G1.
  */

  gsl_matrix *g2covariance;
  /* 
     Covariance matrix for flexion component G2.
  */

  gsl_vector *flexion_redshift;
  /* 
     Flexion redshift vector.
  */

  gsl_vector_int *ccurve;
  /*
    Critical curve estimate vector.
  */

  gsl_vector *ccurve_redshift;
  /* 
     Critical curve estimate redshift vector.
  */
  gsl_vector *ccurveerror;
  /*
    Error on critical curve position.
  */

  gsl_vector *strongfactorBlk;
  /*
    Critical curve prefactor which changes in each inner-level iteration.
  */
  gsl_vector *strongfactorVl;
  /*
    Critical curve prefactor which changes in each inner-level iteration.
    Blk and Vl differ by their power in redshift
  */

  gsl_matrix *msysteminfo;
  /*
    Multiple image system information matrix.
  */


 public:

  MPIDataContainer(ReconstructionOptions &options,FinDifGrid &grid1,int my_rankgiven, int pgiven);
  /*
    Constructor, needs options and a matching finite differencing grid.
  */

  ~MPIDataContainer();
  /*
    Standard destructor, clears (a lot) of memory.
  */

  void readData(ReconstructionOptions &options,FinDifGrid &grid,int iterationindex);
  /*
    Reads the data from FITS files which position are in the options file 
    into process with my_rank 0. This data has to be sent later.
  */
 
  void sendData();
  /*
    Sends the content of the whole DataContainer to all processes.
  */

  void sendResult();
  /*
    Adds together all individual coefficient matrices and data vectors.
    Final result is carried by prcoess with my_rank 0. You have to give
    process rank and total number of processes.
  */
  void resetResult();
  /*
    Set all derived quantities like datavector or coeff to 0 again.
  */

  gsl_vector* showdatavector(string selection);
  /*
    Returns a pointer to the data gsl_vectors. Selection are: ellip1, ellip2,
    shear_redshift, f1, f2, g1, g2, flexion_redshift, ccurve_redshift, 
    ccurveerror, strongfactorBlk, strongfactorVl and datavector.
  */

  gsl_vector_int* showintdatavector(string selection);
  /*
    Returns a pointer to the data gsl_vectors_int. Selection is ccurve.
  */

  gsl_matrix* showdatamatrix(string selection);
  /*
    Returns a pointer to the data gsl_matrices. Selection are: 
    ellip1_covariance, ellip2_covariance, Fij1,Fij2, f1_covariance, 
    f2_covariance, g1_covariance, g2_covariance, msysteminfo, and coeff.
  */
};

#endif 	    /* !MPI_DATACONTAINER_H_ */
