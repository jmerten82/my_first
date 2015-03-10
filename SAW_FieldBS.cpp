#include <iostream>
#include <ctime>
#include <string>
#include <sstream>
#include "options.h"
#include "inputfield.h"
#include "bootstrap.h"

using namespace std;

int main(int argc, char *argv[])
{


  cout <<"----------SAW 1.0 Field BOOTSTRAP routine----------" <<endl;
  cout <<endl;;
  cout <<"Chosen configuration file: " <<argv[1] <<endl;


  BootstrapFieldOptions options1(argv[1]);

  cout <<options1.showflag(1) <<endl;

  gsl_matrix *data;
  string in = options1.showfilename(0);
  ifstream read(in.c_str());
  ifstream input(in.c_str());
  int index = 0;
  double value;
  double value1;
  double value2;
  double value3;
  double value4;
  double value5;


  if(!options1.showflag(1))
    {
      while(read >>value >>value >>value >>value)
	{
	  index++;
	}
      
      read.close();
      data = gsl_matrix_calloc(index,4);
      index = 0;
  while(input >>value1 >>value2 >>value3 >>value4)
    {
      gsl_matrix_set(data,index,0,value1);
      gsl_matrix_set(data,index,1,value2);
      gsl_matrix_set(data,index,2,value3);
      gsl_matrix_set(data,index,3,value4);
      index++;
    }
  input.close();


  Bootstrap bs1(index,4,options1.showbootstrapindex(0),options1.showbootstrapindex(1),options1.showfilename(1)+"/data_realisation");
  bs1.read_data(data);
  bs1.run();
    }

  else
    {
      while(read >>value >>value >>value >>value >>value)
	{
	  index++;
	}
      
      read.close();
      data = gsl_matrix_calloc(index,5);
      index = 0;
      while(input >>value1 >>value2 >>value3 >>value4 >>value5)
    {
      gsl_matrix_set(data,index,0,value1);
      gsl_matrix_set(data,index,1,value2);
      gsl_matrix_set(data,index,2,value3);
      gsl_matrix_set(data,index,3,value4);
      gsl_matrix_set(data,index,4,value5);
      index++;
    }
  input.close();

  Bootstrap bs1(index,5,options1.showbootstrapindex(0),options1.showbootstrapindex(1),options1.showfilename(1)+"/data_realisation");
  bs1.read_data(data);
  bs1.run();
    }


    



  for(int i = options1.showbootstrapindex(0); i <= options1.showbootstrapindex(1); i++)
    {
      ostringstream stream2;
      string in1;
      stream2 << options1.showfilename(1)+"/data_realisation" <<i <<".dat";
      in1 = stream2.str();

      for(int j = 0; j < options1.shownumofsteps(); j++)
	{
	  ostringstream stream3;
	  string out2;
	  stream3 << options1.showfilename(1)+"/ellip_realisation" <<i <<"_" <<options1.showstepbyindex(j)<<".fits";
	  out2 = stream3.str();
	  if(!options1.showflag(1))
	    {

	      cout <<options1.showcuts(7) <<endl;
	      cout <<in1 <<endl;
	      MaskedInputField field1("NA",in1,options1.showstepbyindex(j),options1.showcuts(0),options1.showcuts(1),options1.showcuts(2),options1.showcuts(3),options1.showcuts(4),options1.showcuts(5),options1.showcuts(6),options1.showcuts(7),options1.showflag(0));
	      field1.checkmask(options1.showfilename(2));
	      field1.average(options1.showgal(),options1.showradius());
	      field1.docovariance(options1.showpixelradius());
	      field1.writetofits(out2,0);
	    }
	  else
	    {
	      WeightedInputField field1("NA",in1,options1.showstepbyindex(j),options1.showcuts(0),options1.showcuts(1),options1.showcuts(2),options1.showcuts(3),options1.showcuts(4),options1.showcuts(5),options1.showcuts(6),options1.showcuts(7),options1.showflag(0));
	      field1.checkmask(options1.showfilename(2));
	      field1.average(options1.showgal(),options1.showradius());
	      field1.docovariance(options1.showpixelradius());
	      field1.writetofits(out2,0);
	    }

	  out2 = "";
	}
    }

  return 0;

}




    
