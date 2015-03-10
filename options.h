#ifndef   	OPTIONS_H_
# define   	OPTIONS_H_

#include <iostream>
#include <stdexcept>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <gsl/gsl_matrix.h>
#include <ctime>


using namespace std;

//auxilliary functions

string parametercut(string);  
// cuts out substring between two kinds of given string

int countsymbol(string, const char&); 
// counts the given char in a string


void cutints(string, string, gsl_vector_int*,int);
//gives all integers separated by string in a string to gsl_vector_int 
						
void cutdoubles(string, string, gsl_vector*,int);
//same for doubles

class FieldOptions

/**
   Main options class for the field preparation driver. Once read from
   a parameter file, the following routines read this class for the 
   field preparation options and parameters.
**/

{
 private:

  bool shear;
  /*
    Flag if shear should be used.
  */

  bool flexion;
  /*
    Flag if flexion should be used.
  */

  bool msystems;
  /* 
     Flag if multiple image systems should be used.
  */
  bool lowresccurve;
  /*
    Flag if critical curve estimates should be used.
  */
  
  bool square;
  /*
    Flag if the field should be treated as square with respect to the 
    input x-dimension
  */

  string ellipinput;
  /* 
     Filename of ellipticity catalog
  */

  string flexioninput;
  /*
    Filename of flexion input
  */

  string ccurveinput;
  /*
    Filename of CCurve catalog
  */

  string msysteminput;
  /*
    Filename for multiple image system catalog
  */

  string maskinput;
  /*
    Filename for the mask-defining input
  */

  string outputprefix;
  /*
    Prefix for shear output files
  */

  string flexionprefix;
  /*
    Prefix for the flexion output
   */

  string ccurveprefix;
  /*
    Prefix for the ccurve output
  */

  string msystemprefix;
  /* 
     Prefix for msystem output file
  */

  string clustername;
  /*
    Name for the cluster which appears in the output FITS headers
  */

  gsl_vector *cuts;
  /* 
     Points where the output field should be cut, x1,x2,y1,y2
  */

  gsl_vector *cuts2;
  /* 
     Points where the averaging area should be cut, x1,x2,y1,y2
  */

  gsl_vector_int *steps;
  /* 
     The different resolutions given by start/stopdim and stepsize
  */

  int numgal;
  /*
    Number of galaxies to average over for ellip
  */

  int flexionnumgal;
  /*
    Number of galaxies to average over for flexion
  */

  int startdim;
  /* 
     Lowest field preparation resolution in x-coordinate
  */

  int stopdim;
  /* 
     Highest field preparation resolution in x-coordinate
  */

  int stepsize;
  /*
    stepsize between consecutive resolutions
  */

  double radiusinc;
  /* 
     Radius increment in adaptive averaging routine
  */

  int covpixelradius;
  /*
    Maximum distance to a pixel for which covariances are calculated,
    this mainly defines the speed of the routine but needs good knowledge
    of covariances in the field. Agressive values would be like 2-5,
    conservative like 10-20.
  */

  int nummsystems;
  /* 
     Number of multiple image systems
  */

  int numbootstraps;
  /*
    The number of bootstrap realisations.
  */

  int bootstrapiterator;
  /*
    Shows the position in a bootstrap sequence.
  */

 public:

  FieldOptions(const string& filename);
  /* 
     Constructor, reads all options from filename
  */

  ~FieldOptions();
  /* 
     Standard Destructor, clears (a lot of) memory
  */

  void setbool(string selection, bool value);
  /* 
     Sets one of the boolean values, selections are: shear, flexion, 
     ccurve, msystems, square
  */

  bool showbool(string selection);
  /*
    Returns one of the bools, selections are: shear, flexion, 
    ccurve, msystems, square
  */

  void setstring(string selection, string value);
  /*
    Sets one of the integer values, selections are: ellipinput, flexioninput, 
    ccurveinput, msysteminput, maskinput, shearoutput, flexionoutput, 
    ccurveoutput, msystemoutput, clustername
  */

  string showstring(string selection);
  /*
    Return one of the strings, selections are:
    ellipinput, flexioninput, ccurveinput, msysteminput, 
    maskinput, shearoutput, flexionoutput,ccurveoutput, 
    msystemoutput, clustername
  */ 

  string showstring(string selection,int iterationindex);
  /*
    Returns a full filename not just a prefix depending on the iterationindex.
    Selections are:     ellipinput, flexioninput, ccurveinput, msysteminput, 
    maskinput, shearoutput, flexionoutput,ccurveoutput, 
    msystemoutput, clustername.
    Some of them don't depend then on the index
  */ 

  void setdouble(string selection, double value);
  /*
    Sets one of the double values, selections are: areax1, areax2, areay1, 
    areay2, fieldx1, fieldx2, fieldy1, fieldy2, radiusincrement
  */

  double showdouble(string selection);
  /*
    Returns one of the double values, selections are: areax1, areax2, areay1, 
    areay2, fieldx1, fieldx2, fieldy1, fieldy2, radiusincrement
  */

  void setint(string selection, int value);
  /*
    Sets one of the integer values, selections are: sheargal, flexiongal, 
    startdim, stopdim, stepsize, pixelradius, nummsystems and bootstrap
  */

  int showint(string selection);
  /*
    Returns one of the integer values, selections are: sheargal, flexiongal, 
    startdim, stopdim, stepsize, pixelradius, nummsystems and bootstrap,
    BSiterator
  */

  int showstep(int iterationindex);
  /*
    Returns the x-dim resolution at a given iterationindex of the
    field preparation process.
  */

  int numres();
  /*
    Returns the total number of resolutions which should be processed.
  */

  void printall();
  /*
    Returns  all options to console, mostly for testing
  */

  void iterate();
  /*
    Iterates the bootstrap iterator.
  */

  string show_filename(string selection, int olevel);
  /*
    Returns the right filenames, depending on the bootstrap iterator and 
    the outer level iteration stage. Selections are: shear, flexion, ccurve,
    msystems, ellipinput, ccurveinput, flexioninput, msystemsinput.
  */

};

class ReconstructionOptions
{

  /**
     Defines all the options which are necessary during a reconstruction run.
     This class is used by many reconstruction subroutines since it contains
     most of the parameters used by these routines.
  **/

 private:

  bool shear;
  /*
    Flag indicating if shear should be used in reconstruction
  */

  bool flexion;
  /*
    Flag indicating if flexion should be used in reconstruction
  */

  bool ccurve;
  /*
    Flag indicating if critical curve estimators 
    should be used in reconstruction
  */

  bool msystems;
  /*
    Flag indicating if multiple image systems 
    should be used in reconstruction
  */

  bool weaklimit;
  /*
    Flag indicating if you want to use the assumption that \epsilon = \gamma
    instead of \epsilon = g.
  */

  bool manualsteps;
  /*
    Flag indicating if the reconstruction steps are set manually by the user
  */

  int nummsystems;
  /*
    Number of multiple image systems
  */

  gsl_vector_int *steps;
  /*
    The different reconstruction resolutions stored in a vector
  */

  int startdim;
  /*
    The starting resolution of the reconstruction
  */

  int stopdim;
  /*
    The final resolution in the reconstruction
  */

  int stepsize;
  /*
    The stepsize between different reconstruction resolutions
  */

  int iterations;
  /*
    The maximal number of inner level iterations
  */

  double threshold;
  /*
    The minimal convergence criterium for the inner level iterations
  */

  string shearinput;
  /*
    Filename for the input shear field
  */

  string flexioninput;
  /*
    Filename for the flexion input field
  */

  string ccurveinput;
  /*
    Filename for the critical curve estimator input
  */

  string msysteminput;
  /*
    Filename for the multiple image system input
  */

  string output;
  /*
    Filename for the output files
  */

  string regtype;
  /*
    The type of regularisation used in the reconstruction.
  */

  double regshear;
  /*
    Regularisation parameter for the shear
  */

  double regflexion;
  /*
    Regularisation parameter for flexion
  */

  double shearz;
  /*
    Average or median redshift of the shear background sources
  */

  double flexionz;
  /*
    Avergage or median redshift of the flexion background sources
  */

  double clusterz;
  /*
    Redshift of the reconstructed object
  */

  int numbootstraps;
  /*
    The number of bootstrap realisations.
  */

  int bootstrapiterator;
  /*
    The position in the bootstrap sequence.
  */

 public:

  ReconstructionOptions(const string& filename);
  /*
   Constructor which defines the options from a given input ASCII file filename
  */

  ~ReconstructionOptions();
  /*
    Standard deconstructor, just frees the steps vector and empties the
    filenames
  */

  void setbool(string selection, bool value);
  /*
    Sets one of the flags to value, selections are: shear, weaklimit,
    flexion, ccurve, msystems, manualsteps
  */

  bool showbool(string selection);
  /* 
     Prints a boolean value, selections are: shear, weaklimit,
     flexion, ccurve, msystems, manualsteps
  */

  void setdouble(string selection, double value);
  /* Sets one of the doubles to value, selections are: threshold, 
     regshear, regflexion, shearz, flexionz, clusterz
  */

  double showdouble(string selection);
  /* 
     Prints a double, selection are: threshold, regshear, regflexion, 
     shearz, flexionz, clusterz
  */

  void setint(string selection, int value);
  /* 
     Sets one of the integers to value, selections are: nummsystems, 
     startdim, stopdim, stepsize, iterations, bootstrap
  */

  int showint(string selection);
  /* Prints an integer, selections are: nummsystems, 
     startdim, stopdim, stepsize, iterations, bootstrap, BSiterator
  */

  void setstring(string selection, string value);
  /*
    Sets one of the strings to value, selections are: shearinput, 
    flexioninput, ccurveinput, msystmeinput, output, and regtype
  */

  string showstring(string selection);
  /* 
     Prints a string, selections are: shearinput, flexioninput, 
     ccurveinput, msystmeinput, output and regtype
  */
 

  string showstring(string selection, int iterationindex);
  /*
    Prints a full filename and not just a prefix, according to an 
    iteration index. Selections are: shearinput, flexioninput, ccurveinput,
    msysteminput and output
  */


  void printall();
  /* 
     Prints all options to console out, mainly for testing
  */

  int showstep(int iterationindex);
  /* 
     Shows the reconstruction step vector at a given iterationindex
  */

  int shownumsteps();
  /*
    Return the total number of iteration steps.
  */

  void iterate();
  /*
    Forwards the bootstrapiterator by one.
  */

  string show_filename(string selection, int olevel);
  /*
    Returns filenames depending on the bootstrap sequence positons and
    outer-level iteration sequence. Selections are: shearinput, flexioninput,
    ccurveinput, msystemsinput, output.
  */

};


class HighresOptions
{

 private:

  bool ccurve;
  /*
    Flag indicating if critical curve estimate constraints are used in 
    reconstruction
  */
  bool msystems;
  /*
    Flag indicating if multiple image system constraints are used in 
    reconstruction
  */
  string lowresresult;
  /* 
     Path to the file of the low-resolution result as FITS file.
  */
  string highresccurve;
  /*
    Path to the FITS file containing the critical curve estimates on the 
    final, high resolution.
  */
  string highresmsystemsfile;
  /*
    Path to the ASCII file containing the  multiple image systems on the 
    final, high resolution.
  */
  string interpolationoutput;
  /* 
     Path and filename to the output of the interpolation of the
     low-resolutioon result. Mainly for cross-checking.
  */
  string msystemsoutput;
  /* 
     Path and filename to the ASCII file containing the multiple image
     systems on the cutted cluster core field.
  */ 
  string mergedoutput;
  /*
    Path and filename to the FITS file containing the high-resolution result 
    merged with the low-resolution result.
  */
  int finalx;
  /*
    Resolution in x-direction of the final reconstruction step.
  */
  int finaly;
  /*
    Resolution in y-direction of the final reconstruction step.
  */
  gsl_vector_int *cuts;
  /* 
     Containing x_min,x_max,y_min and y_max to cut the reconstruction field.
  */
  int nummsystems;
  /*
    Number of mupltiple image systems, do not set to 0. Switch off msystems
    with flag above.
  */
  double reg;
  /*
    Regularisation parameter.
  */
  string regscheme;
  /*
    Scheme of regularisation, selections are convergence, convshear, flexion
    and convshearflex.
  */
  double clredshift;
  /*
    Redshift of the lens..
  */
  double ssigma;
  /* 
     Error estimate on the critical curve position.
  */
  int numbootstraps;
  /*
    The number of bootstrap realisations.
  */
  int bootstrapiterator;
  /*
    The position in the bootstrap sequence.
  */

 public:

  HighresOptions(const string&);
  /*
    Constructor, needs position of the parameter file.
  */
  ~HighresOptions();
  /*
    Standard destructor.
  */
  void setbool(bool value,string selection);
  /*
    Sets a flag to value. Selections are ccurve and msystems.
  */ 
  void setstring(string value, string selection);
  /* 
     Sets a string to value, selections are: lowresresult,
     highresccurve, highresmsystemsfile, interpolationoutput, msystemsoutput,
     mergedoutput,regscheme
  */
  void setint(int value, string selection);
  /*
    Sets an integer to value. Selections are: finalx, finaly, x_min,x_max,
    y_min, y_max, bootstrap.
  */
  void setdouble(double value, string selection);
  /*
    Sets a double to value. Selections are: reg, clredshift,ssigma.
  */
  bool showbool(string selection);
  /*
    Prints a boolean value. Selections are: ccurve and msystems.
  */
  string showstring(string selection);
  /*
    Prints a string. Selections are: lowresresult,
    highresccurve, highresmsystemsfile, interpolationoutput, msystemsoutput,
    mergedoutput,regscheme.
  */
  int showint(string selection);
  /*
    Prints an integer value. Selections are: finalx, finaly, x_min,x_max,
    y_min, y_max, bootstrap, BSiterator.
  */
    
  double showdouble(string selection);
  /*
    Prints a double value. Selections are: reg, clredshift,ssigma.
  */

  string show_filename(string filename);
  /*
    Shows the highres in and out filenames, according to the position in the
    bootstrap sequence. Selections are: lowres_in, ccurve_in, msystem_in, 
    msystem_reset, core_out and merge_out.
  */

  void iterate();
  /*
    Iterates the bootstrap sequence.
  */

};

class BootstrapFieldOptions
{

 private:

  int bsstart;
  int bsstop;
  int startres;
  int stopres;
  int stepsize;
  gsl_vector_int *steps;
  bool square;
  bool weight;
  string ellipinput;
  string maskinput;
  string outputdir;
  gsl_vector *cuts1;
  gsl_vector *cuts2;
  int numgal;
  double radiusinc;
  int pixelradius;
  int numofsteps;

 public:

  BootstrapFieldOptions(const string&);
  ~BootstrapFieldOptions();
  int showbootstrapindex(int);
  // 0 for start, 1 for stop
  int showresolution(int);
  // 0 for start, 1 for stop
  int showstepsize();
  int showstepbyindex(int);
  int shownumofsteps();
  bool showflag(int);
  // 0 for square field, 1 for weighted average
  string showfilename(int);
  // 0 for input, 1 for output, 2 for mask
  double showcuts(int);
  // 0,1,2,3 for x1,x2,y1,y2 of area
  // 4,5,6,7 for output field
  int showgal();
  double showradius();
  int showpixelradius();

}; 

class BootstrapReconstructionOptions
{

 private:

  int bsstart;
  int bsstop;
  int startdim;
  int stopdim;
  int stepsize;
  int numofsteps;
  string infile;
  string outfile;
  double regparameter;
  double galredshift;
  double clusterredshift;
  gsl_vector_int *steps;

 public:

  BootstrapReconstructionOptions(const string&);
  ~BootstrapReconstructionOptions();
  int showbootstrapindex(int);
  // 0 for start, 1 for stop
  int showresolution(int);
  // 0 for start, 1 for stop
  int showstepsize();
  int showstepbyindex(int);
  int shownumofsteps();
  string showfilename(int);
  // 0 for input, 1 for output, 2 for mask
  double showreg();
  double showredshift(int);

};

class BootstrapAnalysisOptions
{

 private:

  int realisations;
  string indir;
  string outdir;
  int initres;
  int interres;

 public:

  BootstrapAnalysisOptions(const string&);
  ~BootstrapAnalysisOptions();
  int showrealisations();
  string showfilename(int);
  //0 for input, 1 for output
  int showresolution(int);
  // 0 for init, 1 for interpolated

};

#endif 	    /* !OPTIONS_H_ */
