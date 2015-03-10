#ifndef   	RW_FITS_H_
# define   	RW_FITS_H_


#include <CCfits/CCfits>
#include <valarray>
#include <gsl/gsl_matrix.h>

using namespace std;

/**
   Routines which handle the complete read/write tasks for FITS-files.
**/

void write_pimg(string filename, gsl_matrix *data);
/* 
   Writes a double-valued gsl_matrix data to the primary image of a FITS file
   called filename.
*/

void write_pimgint(string filename, gsl_matrix_int *data);
/* 
   Writes a integer-valued gsl_matrix data to the primary image of a FITS file
   called filename.
*/

void write_pimg(string filename, int x_dim, int y_dim, gsl_vector *data);
/* 
   Writes a double-valued gsl_vector data to the primary image of a FITS file
   called filename. The dimensions of this image are defined by x_dim and
   y_dim, so the vector should have length  (x_dim*y_dim).
*/

void write_pimgint(string filename, int x_dim, int y_dim, gsl_vector_int *data);
/* 
   Writes a integer-valued gsl_vector data to the primary image of a FITS file
   called filename. The dimensions of this image are defined by x_dim and
   y_dim, so the vector should have length  (x_dim*y_dim).
*/

void write_imge(string filename, string extname, gsl_matrix *data);
/* 
   Writes a double-valued gsl_matrix data to the image extension extname
   of a FITS file called filename.
*/

void write_imgeint(string filename, string extname, gsl_matrix_int *data);
/* 
   Writes a integer-valued gsl_matrix data to the image extension extname
   of a FITS file called filename.
*/

void write_imge(string filename, string extname, int x_dim, int y_dim, gsl_vector *data);
/* 
   Writes a double-valued gsl_vector data to an image extension extname 
   of a FITS file called filename. The dimensions of this image are defined 
   by x_dim and y_dim, so the vector should have length  (x_dim*y_dim).
*/

void write_imgeint(string filename, string extname, int x_dim, int y_dim, gsl_vector_int *data);
/* 
   Writes a integer-valued gsl_vector data to an image extension extname 
   of a FITS file called filename. The dimensions of this image are defined 
   by x_dim and y_dim, so the vector should have length  (x_dim*y_dim).
*/

void read_pimg(string filename, gsl_matrix *data);
/*
  Reads a double-valued primary image of a FITS file filename into an 
  extenal gsl_matrix data.
*/

void read_pimgint(string filename, gsl_matrix_int *data);
/*
  Reads an integer-valued primary image of a FITS file filename into an 
  extenal gsl_matrix data.
*/

void read_pimg(string filename, gsl_vector *data);
/*
  Reads a double-valued primary image of a FITS file filename into an 
  extenal gsl_vector data.
*/

void read_pimgint(string filename, gsl_vector_int *data);
/*
  Reads a integer-valued primary image of a FITS file filename into an 
  extenal gsl_vector data.
*/

void read_imge(string filename, string extname, gsl_matrix *data);
/*
  Reads a double-valued image extension extname of a FITS file filename
  into an extenal gsl_matrix data.
*/

void read_imgeint(string filename, string extname, gsl_matrix_int *data);
/*
  Reads a integer-valued image extension extname of a FITS file filename
  into an extenal gsl_matrix data.
*/

void read_imge(string filename, string extname, gsl_vector *data);
/*
  Reads a double-valued image extension extname of a FITS file filename
  into an extenal gsl_vector data.
*/
void read_imgeint(string filename, string extname, gsl_vector_int *data);
/*
  Reads a integer-valued image extension extname of a FITS file filename
  into an extenal gsl_vector data.
*/

void write_header(string filename, string name, double value, string description);
/*
  Writes a double-valued header with its name to the pHDU of a FITS file 
  filename with a description.
*/

void write_header(string filename, string name, int value, string description);
/*
  Writes a integer-valued header with its name to the pHDU of a FITS file 
  filename with a description.
*/
void write_header(string filename, string name, string value, string description);
/*
  Writes a string-valued header with its name to the pHDU of a FITS file 
  filename with a description.
*/

double read_doubleheader(string filename,string name);
/* 
   Reads a double-valued header name from the pHDU of FITS file filename.
*/

int read_intheader(string filename,string name);
/* 
   Reads a integer-valued header name from the pHDU of FITS file filename.
*/

string read_stringheader(string filename,string name);
/* 
   Reads a string-valued header name from the pHDU of FITS file filename.
*/




#endif 	    /* !RW_FITS_H_ */
