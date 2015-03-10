#ifndef   	RECONSTRUCTION_KERNEL_H_
# define   	RECONSTRUCTION_KERNEL_H_

#include <iostream>
#include <gsl/gsl_matrix.h>
#include "galaxycluster.h"
#include "masked_fin_dif.h"
#include "options.h"
#include "util.h"

using namespace std;

/**
   This header collects the routines which are directly involved in the
   reconstruction process. This is the non-MPI version of the reconstruction,
   some MPI versions of these routines can be found in "mpi_coreroutines.h". 
 **/

void StrongFactor(GalaxyCluster &cluster,gsl_vector*redshift,double sigma, gsl_vector *strongfactorBlk,gsl_vector *strongfactorVl);
/*
  Calculates the prefactor matrix for the critical curve term as
  in arXiv:0806.1967. The result is written to strongfactor.
*/

void HighresBlk(HighresOptions &options1,double frct,gsl_vector_int *ccurve, gsl_matrix *msysteminfo,gsl_vector *strongfactor,gsl_matrix *Blk);

void HighresVl(HighresOptions &options1,GalaxyCluster &reference,double frct,gsl_vector_int *ccurve, gsl_matrix *msysteminfo,gsl_vector *strongfactor,gsl_vector *Vl);

#endif 	    /* !RECONSTRUCTION_KERNEL_H_ */
