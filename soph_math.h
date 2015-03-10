#ifndef   	SOPH_MATH_H_
# define   	SOPH_MATH_H_

#include <cmath>
#include "cblas.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "clapack.h"

using namespace std;

/**
   A wrapper to be able to feed gsl matrices into the CLAPACK routines.
**/

void invert_gsl(gsl_matrix *inmatrix);
/*
  Inverts a square gsl_matrix inmatrix and overwrites the original input.
*/

void solve_gsl(gsl_matrix *coeff, gsl_vector *data);
/*
  Solves a linear system of equations with coefficient matrix coeff and result
  vector data by overwriting data
*/



#endif 	    /* !SOPH_MATH_H_ */
