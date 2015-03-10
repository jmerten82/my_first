#ifndef   	BOOTSTRAP_H_
# define   	BOOTSTRAP_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include "rw_fits.h"
#include "interpol.h"
#include "util.h"
#include "cat_writer.h"


using namespace std;

class Bootstrap
{

 private:

  int length;
  int datalength;
  int realisationstart;
  int realisationstop;
  gsl_matrix *data;
  string filename;

 public:

  Bootstrap(int,int,gsl_matrix*,string);
  Bootstrap(int,int,CatWriter*,string);
  ~Bootstrap();
  void run(int sequencenumber);
  void run_delta(int sequencenumber);
  void set_length(int);
  void set_datalength(int);
  void set_realisations(int,int);
  void set_filename(string);
  int show_length();
  int show_datalength();
  int show_realisations(int);
  //0 for startindex, 1 for stop
  string show_filename();

};
/*
class BootstrapAnalysis
{

 private:

  int length;
  int initres;
  int interres;
  gsl_matrix *datamatrix;
  gsl_vector *resultmap;
  gsl_vector *meanfield;
  int outputselection;
  //0 for mean, 1 for standard deviation, 2 for variance
  string inputfile;

 public:

  BootstrapAnalysis(int,int,int,string);
  ~BootstrapAnalysis();
  void analyse(int);
  void write(gsl_vector*);
  void write(string);

};

*/





     



#endif 	    /* !BOOTSTRAP_H_ */
