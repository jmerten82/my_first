#ifndef   	INTERPOL_H_
# define   	INTERPOL_H_

#include <iostream>
#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_statistics.h>
#include "util.h"

using namespace std;

void cspline_mapinterpol(gsl_matrix* in,gsl_matrix *out);
/*
  Interpolates a input map in of certain dimensions to an output map 
  of certain dimensions.
*/

void cspline_mapinterpol_smooth(gsl_matrix* in,gsl_matrix* out);
/*
  Smoothly interpolating the input map in into output map out by just
  increasing the x-dimension by one step in each iteration.
*/

void cspline_maskedinterpol_row(gsl_matrix* in,gsl_matrix_int* inmap,gsl_matrix*out);
/*
  Interpolating a masked map in with mask inmap into out by first interpolating
  the rows of the map.
*/

void cspline_maskedinterpol_col(gsl_matrix *in,gsl_matrix_int *inmap,gsl_matrix*out);
/*
  Interpolating a masked map in with mask inmap into out by first interpolating
  the columns of the map.
*/

void cspline_morecleverinterpol(gsl_matrix *in, gsl_matrix_int *inmask, gsl_matrix *out);
/*
  Routine which extraordinarily clever interpolates in into out given inmask,
  actually uses cspline_cleverinterpol(...) which is not as clever...but this
  one cheats a little by artificially setting pixel values to values which 
  hopefully do not hurt.
*/

void cspline_smooth_morecleverinterpol(gsl_matrix *in, gsl_matrix_int *inmask, gsl_matrix *out);
/*
  Smooth version of cspline_morecleverinterpol.
*/
void cspline_cleverinterpol(gsl_matrix *in, gsl_matrix_int *inmask, gsl_matrix* out, gsl_matrix_int *outmask);
/*
  See cspline_morecleverinterpol(...) but this one does not set any pixel to
  some value but gives sometimes non-reliable results.
*/

void cspline_smooth_cleverinterpol(gsl_matrix *in, gsl_matrix_int *inmask, gsl_matrix* out,gsl_matrix_int *outmask);
/*
  Smooth version of cspline_cleverinterpol.
*/
 





#endif 	    /* !INTERPOL_H_ */
