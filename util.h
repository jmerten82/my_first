#ifndef   	UTIL_H_
# define   	UTIL_H_

#include <iostream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <time.h>
#include <cmath>
#include <stdexcept>
#include "rw_fits.h"
#include "littlelibastro.h"


using namespace std;

/**
   Basic utilities which might be useful while using SaWLens.
**/ 


#define START(a) start = time(NULL); cout << a << "..." << flush;
#define STOP() cout << "finished (" << difftime(time(NULL),start) << " seconds)" <<endl;
#define CLS (cout <<"\033[2J")
#define LOCATE(z,s) (cout <<"\033["<<z <<';' <<s <<'H')
#define CURSORLEFT(n) (cout <<"\033[" <<n <<'D')
#define CURSORRIGHT(n) (cout <<"\033[" <<n <<'C')
#define CURSORUP(n) (cout <<"\033[" <<n <<'A')
#define CURSORDOWN(n) (cout <<"\033[" <<n <<'B')
/*
  Just some cursor positioning escape routines and timing issues.
*/

const double Pi = acos(-1);
/*
  Everybody needs Pi.
*/

void matvec(gsl_matrix *input, gsl_vector *output);
/*
  Converts a gsl_matrix input of row size y_dim and col size x_dim
  into a vector output of length (x_dim*y_dim) in the sense that all 
  row vector are written one by the other.
*/


void matvec(gsl_matrix_int *input, gsl_vector_int *output);
/*
  Converts a gsl_matrix input of row size y_dim and col size x_dim
  into a vector output of length (x_dim*y_dim) in the sense that all 
  row vector are written one by the other.
*/

void vecmat(gsl_vector *input, gsl_matrix *output);
/*
  Converts a gsl_vector of length (x_dim*y_dim) into a gsl_matrix in the
  sense that arrays of length y_dim are taken out of the vector and written
  under each other in the matrix, x_dim times.
*/
 
void vecmat(gsl_vector_int *input, gsl_matrix_int *output);
/*
  Converts a gsl_vector of length (x_dim*y_dim) into a gsl_matrix in the
  sense that arrays of length y_dim are taken out of the vector and written
  under each other in the matrix, x_dim times.
*/

double change(gsl_vector *one, gsl_vector *two);
/*
  Returns the biggest component of the differecne (one - two)
*/

double change(gsl_matrix *one, gsl_matrix *two);
/*
  Returns the biggest component of the differecne (one - two)
*/

void cut(gsl_vector *invector,int x_dim,int y_dim, int x1, int x2, int y1, int y2, gsl_vector *outvector);
/*
  Returns a cut-out of vector invector of dimensions x_dim*y_dim, with cut 
  points x1, x2 and y1, y2 into outvector of dimension ((x2-x1)*(y2-y1)).
*/
 
void cut(gsl_matrix *inmatrix,int x_dim,int y_dim, int x1, int x2, int y1, int y2, gsl_matrix *outmatrix);
/*
  Returns a cut-out of inmatrix of dimensions x_dim, y_dim, with cut 
  points x1, x2 and y1, y2 into outmatrix  of dimension (x2-x1),(y2-y1).
*/

void cut(gsl_vector_int *invector,int x_dim,int y_dim, int x1, int x2, int y1, int y2, gsl_vector_int *outvector); 
/*
  Returns a cut-out of vector invector of dimensions x_dim*y_dim, with cut 
  points x1, x2 and y1, y2 into outvector of dimension ((x2-x1)*(y2-y1)).
*/

void cut(gsl_matrix_int *inmatrix,int x_dim,int y_dim, int x1, int x2, int y1, int y2, gsl_matrix_int *outmatrix);
/*
  Returns a cut-out of inmatrix of dimensions x_dim, y_dim, with cut 
  points x1, x2 and y1, y2 into outmatrix  of dimension (x2-x1),(y2-y1).
*/

void radialprofile(gsl_matrix *input,int midx,int midy,const char *outfile );
/*
  Returns a radial profile of the matrix input in pixel coordinates calculated
  from the centre given by midx and midy. Result is written as ASCII file to 
  outfile.
*/

void zeronorm(gsl_matrix *convinput);
/*
  Normalises a given convergence map in the sense that the pixel with the
  lowest convergence value is set to zero and the rest is mass-sheet-normalised
  accordingly.
*/

void merge_maps(gsl_matrix *coarsemap, gsl_matrix *finemap, int x1pos, int x2pos, int y1pos, int y2pos, gsl_matrix *outmap);
/*
  Merges coarsemap with higher resolved finemap. Both maps have to belong to
  the same field and x1/2pos, y1/2pos give the insertion coordinates in the
  refined output coordinates
*/

double cosmicweight(double lensredshift, double sourceredshift);
/* 
   Returns the cosmic weight function for a sourceredshift and a lensredshift
   as defined in Bartelmann&Schneider 2001.
*/

void ReadMsystemInfo(string input, gsl_matrix *output);
/* 
   Reads the multiple image system information written by the field 
   preparation routines into the msysteminfo matrix used by the reconstruction
   routines.
*/

void WriteMsystemInfo(gsl_matrix *input,string output);
/*
  Writes a gsl_matrix into a ASCII file.
*/

void read_doubles(string input, int number, vector<double>&);
/*
  Reads a number of double numbers from a string. Subsequent numbers have
  to be separate by any non-number symbol exlcuding . and -.
*/

void read_ints(string input, int number, vector<int>&);
/*
  Reads a number of integer numbers from a string. Subsequent numbers have
  to be separate by any non-number symbol exlcuding . and -.
*/


void read_doubles(string input, vector<double>&);
/*
  Reads all double numbers from a string. Subsequent numbers have
  to be separate by any non-number symbol exlcuding . and -.
*/

void read_ints(string input, vector<int>&);
/*
  Reads all integer numbers from a string. Subsequent numbers have
  to be separate by any non-number symbol exlcuding . and -.
*/

double read_double(string input);
/*
  Reads the first double value from a string and returns it.
*/

int read_int(string input);
/*
  Reads the first integer value from a string and returns it.
*/

string read_word(string input, string symbol);
/*
  Cuts a string between certain markers symbol and returns it.
*/

bool read_mind(string input);
/*
  Try to interpret answer to a question.
*/

void read_sequence(string input, vector<int>&);
/*
  Reads an integer sequence of the form: "a-b,c,d,e-f" and writes
  it into a vector with the full sequence length.
*/


#endif 	    /* !UTIL_H_ */
