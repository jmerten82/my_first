//FINAL MPI
//soph_math.h
//Contains high performance matrix inversion routine based on lapack
//used with ATLAS



#include "soph_math.h"

using namespace std;

void invert_gsl(gsl_matrix *inmatrix)
{

  int *piv = new int[inmatrix->size1];
  clapack_dgetrf(CblasRowMajor,inmatrix->size1,inmatrix->size1,gsl_matrix_ptr(inmatrix,0,0),inmatrix->size1,piv);
  clapack_dgetri(CblasRowMajor,inmatrix->size1,gsl_matrix_ptr(inmatrix,0,0),inmatrix->size1,piv);

}

void solve_gsl(gsl_matrix *coeff, gsl_vector *data)
{

  int *piv = new int[data->size];
  clapack_dgetrf(CblasRowMajor,data->size,data->size,gsl_matrix_ptr(coeff,0,0),data->size,piv);
  clapack_dgetrs(CblasRowMajor,CblasNoTrans,data->size,1,gsl_matrix_ptr(coeff,0,0),data->size,piv,gsl_vector_ptr(data,0),data->size);

}
  


