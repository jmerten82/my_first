/***
    NON-public header file.
    Author: Julian Merten
    Inst.: ITA/ZAH Heidelberg
    2010
***/



#ifndef   	CAT_WRITER_H_
# define   	CAT_WRITER_H_


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <cat_reader.h>


using namespace std;

class CatWriter
{

  /**
    A class which takes an ASCII catalogue read by another cat_reader class
    and produces an output of that.
  **/

  CatReader *input;
  /*
    The underlying cat_reader class which holds data that can be processed
    by the cat_writer.
  */

  vector<bool> row_flags;
  /*
    A vector with the size of the row-size of the catalogue flagging 
    switched-on and turned-off rows.
  */

  vector<int> col_flags;
  /*
    A vector with the size of the column-size of the catalogue flagging 
    switched-on and turned-off columns.
  */

  int num_rows;
  /*
    The number of output rows.
  */

  int num_columns;
  /*
    The number of output columns.
  */

  string sep;
  /*
    The separation symbol in the ASCII file separating the data values.
  */

  string mode;
  /*
    Defines the output mode of the class. Selection are "normal",
    "scientific" and "int".
  */

  bool write_info;
  /*
    Flag determining if info about the columns is written to the ASCII output.
  */

  gsl_matrix *data;
  /*
    Matrix containing the data of the catalogue, but only if access to it is
    necesssary. Relates to alloc.
  */

  bool alloc;
  /*
    Simple flag, showing if memory for data was allocated.
  */

  bool update_needed;
  /*
    Flag which keeps track of the status for the content of data.
  */

  vector<std::string> names;
  /*
    Holds the columns names of the writer. Somewhat ugly construction still.
  */

  bool WCS;
  /*
    Flag, showing if the catalogue contains objects in a WCS frame.
  */


 public:

  CatWriter(CatReader*);
  /*
    Standard constructor which just needs a cat_reader as input. Some
    checks on the cat_reader are performed. With this constructor all
    columns and rows of the cat_reader are activated.
  */

  CatWriter(CatReader*, string col_select);
  /*
    Standard constructor which just needs a cat_reader as input. Some
    checks on the cat_reader are performed. With this constructor only those
    columns of the cat_reader are activated, which fit the selection
    scheme. e.g.: 1-5,7,6,8-25. All rows are activated.
  */

  CatWriter(CatReader*, string col_select, string row_select);
  /*
    Standard constructor which just needs a cat_reader as input. Some
    checks on the cat_reader are performed. With this constructor only those
    columns of the cat_reader are activated, which fit the selection
    scheme. e.g.: 1-5,7,6,8-25. The same applies to the rows. col_select
    can be set to "all".
  */

  ~CatWriter();
  /*
    Standard destructor, freeing memory.
  */

  void set_col(int,int);
  /*
    Turns of or switches on a certain column for further processing.
  */

  void set_col(string);
  /*
    Sets all the output columns due to a certain scheme. e.g. "1-4,6,5,7-100"
  */

  int show_col(int);
  /*
    Shows if a certain column is switched on or turned off and to what value.
  */

  void set_row(int,bool);
  /*
    Turns of or switches on a certain row for further processing.
  */

  void set_row(string);
  /*
    Sets all the output rows due to a certain scheme. e.g. "1-4,6,5,7-100"
  */

  bool show_row(int);
  /*
    Shows if a certain row is switched on or turned off.
  */

  int show_rowdim();
  /*
    Returns the number of rows of the CatWriter.
  */

  int show_coldim();
  /*
    Returns the number of columns of the CatWriter.
  */


  void set_mode(string symbol, string selection);
  /*
    Sets the output mode of the ASCII routine. The symbol defines the
    separation between values. Selections are: normal and scientific. 
  */

  void set_infomode(bool);
  /*
    Sets the flag if an additional info line about the column 
    names is written to the ASCII out.
  */

  string show_mode();
  /*
    Shows the output mode of the ASCII routine.
  */

  bool show_infomode();
  /*
    Shows if the infomode is toggled.
  */

  void write(const string& filename);
  /*
    Writes the data to ASCII.
  */

  void write(gsl_matrix *output);
  /*
    Writes the data into external gsl matrix.
  */

  void create_copy();
  /*
    Creates a physical copy of the actual data content of the cat_writer
    from the cat_reader and writes it into data.
  */

  gsl_matrix* show_data();
  /*
    Returns a pointer to the data content.
  */

  void set_WCS(bool);
  /*
    Turns the WCS system on or off.
  */

  bool show_WCS();
  /*
    True if WCS is activated, false if not.
  */

  void write_DS9_points(const string &filename, int col1, int col2, string selection);
  /*
    Writes a DS9 region file with the shape selection, col1 will be the inter
    preted as x-coordinate and col2 as y-coordinate. Selections are:
    circle, box, diamond, cross, xpoint or boxcircle.
  */

  void write_DS9_circles(const string &filename, int col1, int col2, int col3, double size);
    /*
      Writes a DS9 region file. Objects are circles of radius size, with the
      ids given by col1 and x,y given by col1 and col2.
    */ 

  void add_tocol(int,double);
  /*
    Adds a double value to a certain col.
  */

  void multiply_tocol(int,double);
  /*
    Muliplies a certain column with a double value.
  */

  void power_tocol(int,double);
  /*
    Takes a certain column to a certain power.
  */

  void trig_tocol(int,string selection);
  /*
    Applies a trigonometric function to a certain column. Selections are:
    sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, exp.
  */

  void addfunccol_tocol(int, int);
  /*
    Adds  a certain column to another column.
  */

  void addfunccol_tocol(int, int,double);
  /*
    Adds the power of a certain column to another column.
  */

  void addfunccol_tocol(int, int,string);
  /*
    Adds the function of a certain column to another column. Selections are:
    sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, exp.   
  */

  void multiplyfunccol_tocol(int, int);
  /*
    Multiplies a certain column with another column.
  */

  void multiplyfunccol_tocol(int, int,double);
  /*
    Multiplies the power of a certain column with another column.
  */

  void multiplyfunccol_tocol(int, int,string);
  /*
    Multiplies the function of a certain column to another column. 
    Selections are:
    sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, exp.   
  */

  void kill_row(int col, double thresh1, double thresh2);
  /*
    Sorts out all rows where the values of a certain column do not lie
    within a window defined by thresh1 and thresh2. Be aware that this will not
    work on the original catalogue, but on the modified one.
  */

  void add_col(int pos, double value,string name);
  /*
    Adds a uniform column with a certain value. A name of the column has to be 
    given.
  */  

};


#endif 	    /* !CAT_WRITER_H_ */
