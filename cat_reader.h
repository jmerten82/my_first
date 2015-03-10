/***
    NON-public header file.
    Author: Julian Merten
    Inst.: ITA/ZAH Heidelberg
    2010
***/



#ifndef   	CAT_READER_H_
# define   	CAT_READER_H_


#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <gsl/gsl_matrix.h>
#include <util.h>


using namespace std;

class CatReader
{

  /**
     This class reads an parts of an ASCII file into a GSL matrix, while
     getting rid of comments and unwanted lines. 
  **/

 protected:


  gsl_matrix *data;
  /*
    The matrix containing the actual data.
  */

  vector<std::string> col_names;
  /*
    Contains the names of the different columns in the cat_reader. The zero'th 
    element is the name of the cat_reader.
  */

  int cat_counter;
  /*
    Countint of how many individual catalogues the cat-reader consists.
  */


 public:

  CatReader(const string &filename);
  /*
    Constructor, automatically reading an input ASCII file without
    any limits.
  */

  CatReader(const string &filename, string symbol);
  /*
    Constructor, automatically reading an input ASCII file without
    any limits. Also tries to find the column names, which are indicated 
    in the input by symbol.
  */

  CatReader(const string &filename, int limit1, int limit2);
  /*
    Constructor, automatically reading an input ASCII file. limit1 limits
    the values of read-in lines and limit2 the number of read in columns.
  */

  CatReader(const string &filename, int limit1, int limit2, string symbol);
  /*
    Constructor, automatically reading an input ASCII file. limit1 limits
    the values of read-in lines and limit2 the number of read in columns.
    Also tries to find the column names, which are indicated 
    in the input by symbol.
  */

  ~CatReader();
  /*
    Standard destructor.
  */

  int show_rowdim();
  /*
    Returns the number rows in the read catalogue.
  */

  int show_coldim();
  /*
    Returns the number of columns in the read catalogue.
  */

  void read(const string &filename);
  /*
    Reads another input catalogue into the data of the class. Checks if the
    catalogue row-sizes match are performed.
  */

  void read(const string &filename,string symbol);
  /*
    Reads another input catalogue into the data of the class. Checks if the
    catalogue row-sizes match are performed. Also tries to find the column
    names, indicated in the input by symbol.
  */

  void read(const string &filename, int limit1, int limit2);
  /*
    Reads another input catalogue into the data of the class. Checks if the
    catalogue row-sizes match are performed.
  */

  void read(const string &filename, int limit1, int limit2, string symbol);
  /*
    Reads another input catalogue into the data of the class. Checks if the
    catalogue row-sizes match are performed. Also tries to find the column
    names, indicated in the input by symbol.
  */

  void set_names(int col, string name);
  /*
    Sets the name of the col'th column of the cat_reader. The zero'th column
    is the name of the catalogue.
  */

  string show_names(int col);
  /*
    Shows the name of the col'th column of the cat_reader. The zero'th column
    is the name of the cat_reader.
  */


  int show_counter();
  /*
    Returns the catalogue counter.
  */

  gsl_matrix* show_data();
  /*
    Returns a pointer to the data matrix of the class.
  */

  void write_data(gsl_matrix* output);
  /*
    Writes the data matrix into the output.
  */ 

};



#endif 	    /* !CAT_READER_H_ */



