#ifndef   	MASKED_FIN_DIF_H_
# define   	MASKED_FIN_DIF_H_

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <gsl/gsl_matrix.h>
#include "util.h"
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

class FinDifGrid
{
  /** 
      The essential representation of a finite differencing grid, these
      grids can handle derivatives up to the 3rd order and contain extremely
      fast matrix-multiplication routines. The grid itself needs a negligible
      amount of memory and can also be masked.
  **/


 private:

  double x_frct;
  /*
    if the grid is just a cut-out from a bigger field, this value descibes
    the ratio of the respective x-dimensions.
  */

  int x_dim;
  /*
    The x-dimension of the grid.
  */

  int y_dim;
  /*
    The y-dimension of the grid.
  */

  int fieldpixels;
  /*
    The number of non-masked pixels in the grid.
  */

  int badpixels;
  /*
    Number of pixels which had to be refused after grid construction
    because the finite differencing schemes would not fit.
  */

  gsl_matrix_int *gridmask;
  /*
    The mask which defines the grid. Format is value -1 for masked pixels and
    an integer number from 0 to fieldpixels for useable pixels. This number
    is there index in the gridvector.
  */

  gsl_matrix_int *typemap;
  /* 
     Matrix which defines all properties of an pixel. For further documentation
     look in the additional material.
  */

  gsl_vector_int *typevector;
  /* Essentially the same as the typemap but contains only the non-masked pixels
     for internal use.
  */

  gsl_matrix_int *distancetable;
  /*
    Table which contains the distances in the fieldvector between neighbouring
    pixels. This is crucial for masked fields. For further reading see the 
    additional material.
  */

  bool thirdorder;
  /* 
     Flag which is true if the grid needs to be able to do third order 
     derivatives. If this is not needed set it to false and the accuracy of
     second order derivatives will slightly improve at the grid borders.
  */

 public:

  FinDifGrid(gsl_matrix_int* fieldmask, double x_frct, bool thirdorder);
  /*
    Constructor, which needs a fieldmask, where value 1 means pixel is
    masked and everything else means non-masked pixel. 
    The x_frct gives the ratio between x-dimension of the whole field 
    and the grid and if thirorder is  true the grid can perform third order 
    derivatives. Set to false if this is not needed.
  */


  ~FinDifGrid();
  /*
    Standard destructor, deallocates the gsl_vectors and matrices.
  */

  int showint(string selection);
  /*
    returns the integer values of the grid, selections are: 
    "x_dim", "y_dim", "dim (x_dim*y_dim)", "fieldpixels", "badpixels"
  */
  
  double showdouble(string selection);
  /*
    returns the double values of the grid, only selection right now is
    "x_frct"
  */

  void givegridmask(gsl_matrix_int* output);
  /*
    returns the fieldmask of the grid to an external gsl_matrix of
    row size y_dim and col size x_dim.
  */

  void givetypemap(gsl_matrix_int *output);
  /*
    returns the typemap of the grid to an external gsl_matrix of
    row size y_dim and col size x_dim.
  */
  
  void givetypemap(gsl_vector_int *output);
  /*
    returns the typemap of the grid to an external gsl_vector of
    length fieldpixels.
  */

  void givedistancetable(gsl_matrix_int *output);
  /*
    returns the distancetable of the grid to an external gsl_matrix of
    row size fieldpixels and col size 10.
  */

  void givedistancetable(const string &filename);
  /*
    writes the distancetable of the grid into an ASCII file named filename
  */


  void maptogridvector(gsl_matrix *input, gsl_vector *output);
  /*
    converts a map of row size y_dim and col size x_dim into a 
    gridvector of length fieldpixels.
  */

  void maptogridvector(gsl_matrix_int *input, gsl_vector_int *output);
  /*
    converts a map of row size y_dim and col size x_dim into a 
    gridvector of length fieldpixels.
  */

  void gridvectortomap(gsl_vector *input, gsl_matrix *output, double maskvalue);
  /* 
     converts a gridvector of length fieldpixels into a map of row size y_dim 
     and col size x_dim by setting masked pixels to maskvalue.
  */

  void gridvectortomap(gsl_vector_int *input, gsl_matrix_int *output, int maskvalue);
  /* 
     converts a gridvector of length fieldpixels into a map of row size y_dim 
     and col size x_dim by setting masked pixels to maskvalue.
  */

  double a1value(int l ,int k);
  /*
    returns the value of the deflection angle 1 finite-differences
    representation at row l and col k of this matrix.
  */


  double a2value(int l ,int k);
  /*
    returns the value of the deflection angle 2 finite-differences
    representation at row l and col k of this matrix.
  */


  double s1value(int l ,int k);
  /*
    returns the value of the shear 1 finite-differences
    representation at row l and col k of this matrix.
  */
 
  double s2value(int l ,int k);
  /*
    returns the value of the shear 2 finite-differences
    representation at row l and col k of this matrix.
  */

  double cvalue(int l ,int k);
  /*
    returns the value of the convergence finite-differences
    representation at row l and col k of this matrix.
  */

  double f1value(int l ,int k);
  /*
    returns the value of the flexion F1 finite-differences
    representation at row l and col k of this matrix.
  */

  double f2value(int l ,int k);
  /*
    returns the value of the flexion F2 finite-differences
    representation at row l and col k of this matrix.
  */

  double g1value(int l ,int k);
  /*
    returns the value of the flexion G1 finite-differences
    representation at row l and col k of this matrix.
  */

  double g2value(int l ,int k);
  /*
    returns the value of the flexion G2 finite-differences
    representation at row l and col k of this matrix.
  */

  void writea1(gsl_matrix *output);
  /*
    writes the deflextion angle 1 finite-differences
    representation to an external gsl_matrix output of row size
    fieldpixels and col size fieldpixels.
  */

  void writea2(gsl_matrix *output);
  /*
    writes the deflextion angle 2 finite-differences
    representation to an external gsl_matrix output of row size
    fieldpixels and col size fieldpixels.
  */

  void writes1(gsl_matrix *output);
  /*
    writes the shear 1 finite-differences
    representation to an external gsl_matrix output of row size
    fieldpixels and col size fieldpixels.
  */
  
  void writes2(gsl_matrix *output);
  /*
    writes the shear 2 finite-differences
    representation to an external gsl_matrix output of row size
    fieldpixels and col size fieldpixels.
  */

  void writec(gsl_matrix *output);
  /*
    writes the convergence finite-differences
    representation to an external gsl_matrix output of row size
    fieldpixels and col size fieldpixels.
  */

  void writef1(gsl_matrix *output);
  /*
    writes the flexion F1 finite-differences
    representation to an external gsl_matrix output of row size
    fieldpixels and col size fieldpixels.
  */

  void writef2(gsl_matrix *output);
  /*
    writes the flexion F2 finite-differences
    representation to an external gsl_matrix output of row size
    fieldpixels and col size fieldpixels.
  */

  void writeg1(gsl_matrix *output);
  /*
    writes the flexion G1 finite-differences
    representation to an external gsl_matrix output of row size
    fieldpixels and col size fieldpixels.
  */

  void writeg2(gsl_matrix *output);
  /*
    writes the flexion G2 finite-differences
    representation to an external gsl_matrix output of row size
    fieldpixels and col size fieldpixels.
  */

  int firstorderscheme(int row , int index);
  /*
    complicated thing. It is a set of iterationindices which makes sure
    that a first order fin-dif matrix mutliplcication at a given row 
    index is done in an optimal way. 
    This set contains 9 indices individually accessable through index.
  */ 

  int secondorderscheme(int row , int index);
  /*
    complicated thing. It is a set of iterationindices which makes sure
    that a second order fin-dif matrix mutliplcication at a given row 
    index is done in an optimal way. 
    This set contains 13 indices individually accessable through index.
  */

  int thirdorderscheme(int row , int index);
  /*
    complicated thing. It is a set of iterationindices which makes sure
    that a third order fin-dif matrix mutliplcication at a given row 
    index is done in an optimal way. 
    This set contains 25 indices individually accessable through index.
  */

  void a1multvec(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection angle 1 fin-dif matrix with a vector
    of length fieldpixels and returns the result of length fieldpixels.
  */

  void a1multvec_fast(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection angle 1 fin-dif matrix with a vector
    of length fieldpixels in a speed-optiomized way
    and returns the result of length fieldpixels.
  */

  void a2multvec(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection angle 2 fin-dif matrix with a vector
    of length fieldpixels and returns the result of length fieldpixels.
  */

  void a2multvec_fast(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection angle 2 fin-dif matrix with a vector
    of length fieldpixels in a speed-optiomized way
    and returns the result of length fieldpixels.
  */

  void s1multvec(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection shear 1 fin-dif matrix with a vector
    of length fieldpixels and returns the result of length fieldpixels.
  */

  void s1multvec_fast(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection shear1 fin-dif matrix with a vector
    of length fieldpixels in a speed-optiomized way
    and returns the result of length fieldpixels.
  */

  void s2multvec(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection shear 2 fin-dif matrix with a vector
    of length fieldpixels and returns the result of length fieldpixels.
  */

  void s2multvec_fast(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection shear 2 fin-dif matrix with a vector
    of length fieldpixels in a speed-optiomized way
    and returns the result of length fieldpixels.
  */

  void cmultvec(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection convergence fin-dif matrix with a vector
    of length fieldpixels and returns the result of length fieldpixels.
  */

  void cmultvec_fast(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection convergence fin-dif matrix with a vector
    of length fieldpixels in a speed-optiomized way
    and returns the result of length fieldpixels.
  */

  void f1multvec(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection flexion F1 fin-dif matrix with a vector
    of length fieldpixels and returns the result of length fieldpixels.
  */

  void f1multvec_fast(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection flexion F1 fin-dif matrix with a vector
    of length fieldpixels in a speed-optiomized way
    and returns the result of length fieldpixels.
  */

  void f2multvec(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection flexion F2 fin-dif matrix with a vector
    of length fieldpixels and returns the result of length fieldpixels.
  */

  void f2multvec_fast(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection flexion F2 fin-dif matrix with a vector
    of length fieldpixels in a speed-optiomized way
    and returns the result of length fieldpixels.
  */

  void g1multvec(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection flexion G1 fin-dif matrix with a vector
    of length fieldpixels and returns the result of length fieldpixels.
  */

  void g1multvec_fast(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection flexion G1 fin-dif matrix with a vector
    of length fieldpixels in a speed-optiomized way
    and returns the result of length fieldpixels.
  */

  void g2multvec(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection flexion G2 fin-dif matrix with a vector
    of length fieldpixels and returns the result of length fieldpixels.
  */
  void g2multvec_fast(gsl_vector *input, gsl_vector *output);
  /*
    multiplies the deflection flexion G2 fin-dif matrix with a vector
    of length fieldpixels in a speed-optiomized way
    and returns the result of length fieldpixels.
  */

  double ShearBlk(int l, int k, gsl_vector *ellip1, gsl_vector *ellip2, gsl_matrix *F1, gsl_matrix *F2,string selection);
  /*
    returns the value of the shear part of the reconstruction coefficient 
    matrix at row l and col k. As input the measured ellipticities and the
    according prefactor matrices are needed. Selection are "shear" as an
    initial step and "reduced shear" after the initial reconstruction.
  */

  double ShearVl(int l, gsl_vector *ellip1, gsl_vector * ellip2, gsl_matrix *F1, gsl_matrix * F2, gsl_vector *redshift, string selection);
  /*
    returns the value of the shear part of the reconstruction data 
    vector at position l. As input the measured ellipticities, their
    redshifts and the according prefactor matrices are needed. 
    Selection are "shear" as an initial step and "reduced shear" after the 
    initial reconstruction.
  */

  double FlexionBlk(int l, int k,gsl_matrix *F1, gsl_matrix *F2,gsl_matrix *G1, gsl_matrix *G2);
  /*
    returns the value of the flexion part of the reconstruction coefficient 
    matrix at row l and col k. As input the prefactor matrices for the four
    flexion components are needed.
  */

  double FlexionVl(int l, gsl_vector *f1, gsl_vector *f2, gsl_vector *g1, gsl_vector *g2, gsl_matrix *F1, gsl_matrix *F2, gsl_matrix *G1, gsl_matrix *G2, gsl_vector *redshift);
  /*
    returns the value of the flexion part of the reconstruction data 
    vector at position l. As input the measured flexions, their
    redshifts and the according prefactor matrices are needed. 
  */
  

  double CcurveBlk(int l, int k, gsl_vector_int *ccurve, gsl_vector *factors);
  /*
    returns the value of the critical curve part of the reconstruction 
    coefficient matrix at row l and col k. As input the critical curve 
    and the according prefactor matrix is needed.
  */

  double CcurveVl(int l, gsl_vector_int *ccurve,gsl_vector *factors);
  /*
    returns the value of the critical curve part of the reconstruction 
    data vector at pos l. As input the critical curve 
    and the according prefactor matrix is needed.
  */

  double MsystemBlk(int l, int k, gsl_matrix *msysteminfo);
  /*
    returns the value of the multiple image systems part of the reconstruction 
    coefficient matrix at row l and col k. As input the matrix which contains
    all msystems reconstruction info is needed.
  */

  double MsystemVl(int l, gsl_matrix *msysteminfo);
  /*
    returns the value of the multiple image systems part of the reconstruction 
    data vector at pos l. As input the matrix which contains
    all msystems reconstruction info is needed.
  */

  double RegularisationBlk(int l, int k, double eta1, double eta2, string selection);
  /*
    returns the value of the regularisation part of the reconstruction 
    coefficient matrix at row l and col k. As input the regularisation 
    parameters for shear and flexion are needed. The selections are:
    "convergence", "convshear", "flexion", "convshearflex" if you want to
    regularise on everything, convergence shear and flexion. 
  */

  double RegularisationVl(int l, gsl_vector *conv, gsl_vector *shear1, gsl_vector *shear2, gsl_vector *f1, gsl_vector *f2, gsl_vector *g1, gsl_vector *g2, double eta1, double eta2, string selection);
  /*
    returns the value of the regularisation part of the reconstruction 
    data vector at pos l. As input the reference convergence, shear and 
    flexion values and the  accroding regularisation parameters are needed. 
    The selections are: "convergence", "convshear", "flexion" 
    or anything else if you want to regularise on everything convergence, 
    shear and flexion. 
  */


};

#endif 	    /* !MASKED_FIN_DIF_H_ */
