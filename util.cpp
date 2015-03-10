//FINAL MPI
//util.h
//some little utility routines
//some of them are very important



#include "util.h"


using namespace std;


//useful time operations and cursor controls

#define START(a) start = time(NULL); cout << a << "..." << flush;
#define STOP() cout << "finished (" << difftime(time(NULL),start) << " seconds)" <<endl;
#define CLS (cout <<"\033[2J")
#define LOCATE(z,s) (cout <<"\033["<<z <<';' <<s <<'H')
#define CURSORLEFT(n) (cout <<"\033[" <<n <<'D')
#define CURSORRIGHT(n) (cout <<"\033[" <<n <<'C')
#define CURSORUP(n) (cout <<"\033[" <<n <<'A')
#define CURSORDOWN(n) (cout <<"\033[" <<n <<'B')


void matvec(gsl_matrix *input, gsl_vector *output)
{
  int counter = 0;
  for(int i = 0; i < input->size1; i++)
    {
      for(int j = 0; j < input->size2; j++)
	{
	  gsl_vector_set(output,counter,gsl_matrix_get(input,i,j));
	  counter++;
	}
    }
}
void matvec(gsl_matrix_int *input, int x_dim, int y_dim, gsl_vector_int *output)
{
  int counter = 0;
  for(int i = 0; i < input->size1; i++)
    {
      for(int j = 0; j < input->size2; j++)
	{
	  gsl_vector_int_set(output,counter, gsl_matrix_int_get(input,i,j));
	  counter++;
	}
    }
}

void vecmat(gsl_vector *input, gsl_matrix *output)
{

  int counter = 0;
  for(int i = 0; i < output->size1; i++)
    {
      for(int j = 0; j < output->size2; j++)
	{
	  gsl_matrix_set(output, i,j, gsl_vector_get(input,counter));
	  counter++;
	}
    }
}

void vecmat(gsl_vector_int *input, gsl_matrix_int *output)
{

  int counter;
  for(int i = 0; i < output->size2; i++)
    {
      for(int j = 0; j < output->size1; j++)
	{
	  gsl_matrix_int_set(output, i,j, gsl_vector_int_get(input, counter));
	  counter++;
	}
    }
}


double change(gsl_vector *one, gsl_vector *two)
{
  gsl_vector_sub(one,two);
  for(int i = 0; i < one->size; i++)
    {
      gsl_vector_set(one, i, abs(gsl_vector_get(one,i)));
    }
  return gsl_vector_max(one);
}

double change(gsl_matrix *one, gsl_matrix *two)
{
  //same for gsl-matrices
  gsl_matrix_sub(one,two);
  for(int i = 0; i < one->size1; i++)
    {
      for(int j = 0; j < one->size2; j++)
	{
	  gsl_matrix_set(one,i,j,abs(gsl_matrix_get(one,i,j)));
	}
    }
  return gsl_matrix_max(one);
}

void cut(gsl_vector *invector,int x_dim,int y_dim, int x1, int x2, int y1, int y2, gsl_vector *outvector)

     //does a selected cutout of a given gsl-vector
{
  int x = x2 - x1;
  int y = y2 - y1;

  for(int i = 0; i < y; i++)
    {
      for(int j = 0; j < x; j++)
	{
	  gsl_vector_set(outvector, i*x+j,gsl_vector_get(invector, (y1+i)*x_dim+(x1+j)));
	}
    }
}

void cut(gsl_matrix *inmatrix,int x_dim,int y_dim, int x1, int x2, int y1, int y2, gsl_matrix *outmatrix)
{

  //same for matrix
  int x = x2 - x1;
  int y = y2 - y1;

  for(int i = 0; i < y; i++)
    {
      for(int j = 0; j < x; j++)
	{
	  gsl_matrix_set(outmatrix, i,j,gsl_matrix_get(inmatrix, y1+i,x1+j));
	}
    }
}

void cut(gsl_vector_int *invector,int x_dim,int y_dim, int x1, int x2, int y1, int y2, gsl_vector_int *outvector)

     //does a selected cutout of a given gsl-vector
{
  int x = x2 - x1;
  int y = y2 - y1;

  for(int i = 0; i < y; i++)
    {
      for(int j = 0; j < x; j++)
	{
	  gsl_vector_int_set(outvector, i*x+j,gsl_vector_int_get(invector, (y1+i)*x_dim+(x1+j)));
	}
    }
}

void cut(gsl_matrix_int *inmatrix,int x_dim,int y_dim, int x1, int x2, int y1, int y2, gsl_matrix_int *outmatrix)
{

  //same for matrix
  int x = x2 - x1;
  int y = y2 - y1;

  for(int i = 0; i < y; i++)
    {
      for(int j = 0; j < x; j++)
	{
	  gsl_matrix_int_set(outmatrix, i,j,gsl_matrix_int_get(inmatrix, y1+i,x1+j));
	}
    }
}

void radialprofile(gsl_matrix *input,int midx,int midy,const char *outfile )
{

  //does a radial profile of a given input matrix wrt a given center
  //and writes it into ascii-file
  int x_dim = input->size2;
  int y_dim = input->size1;
  int dim = x_dim*y_dim;
  gsl_matrix *workingmatrix = gsl_matrix_calloc(dim,2);
  gsl_matrix *result = gsl_matrix_calloc(dim,2);

  for(int i = 0; i < y_dim; i++)
    {
      for(int j = 0;j < x_dim; j++)
	{
	  gsl_matrix_set(workingmatrix,i*x_dim+j,1,(i-midy)*(i-midy)+(j-midx)*(j-midx));
	  gsl_matrix_set(workingmatrix,i*x_dim+j,0,gsl_matrix_get(input,i,j));
	}
    }

  for(int i = 0;i < dim; i++)
    {
      double value = 0.0;
      int counter = 0;

      for(int j = 0;j < dim; j++)
	{
	  if(gsl_matrix_get(workingmatrix,j,1)==gsl_matrix_get(workingmatrix,i,1))
	    {
	      value = value + gsl_matrix_get(workingmatrix,j,0);
	      counter++;
	    }
	}
      value = value/counter;
      gsl_matrix_set(result,i,0,gsl_matrix_get(workingmatrix,i,1));
      gsl_matrix_set(result,i,1,value);
    }  
  ofstream output(outfile);
  for(int i = 0; i < dim; i++)
    {
      output <<sqrt(gsl_matrix_get(result,i,0)) <<"\t" <<gsl_matrix_get(result,i,1) <<endl;
    }
}

void zeronorm(gsl_matrix *convinput)
{
  //normalises a given conv-map to 0 wrt to faintest point applying
  //mass-sheet-degenerancie-trafo 
  double minimum = gsl_matrix_min(convinput);
  double lambda = 1.0/(1.0-minimum);
  gsl_matrix_scale(convinput,lambda);
  gsl_matrix_add_constant(convinput,1.0-lambda);
}
void zeronorm(gsl_vector *convinput)
{
  double minimum = gsl_vector_min(convinput);
  double lambda = 1.0/(1.0-minimum);
  gsl_vector_scale(convinput,lambda);
  gsl_vector_add_constant(convinput,1.0-lambda);
}

void merge_maps(gsl_matrix *coarsemap, gsl_matrix *finemap,int x1pos, int x2pos, int y1pos, int y2pos, gsl_matrix *outmap)
{

  //merges a high and a low resolution gsl map, you need to give both
  //maps and their resolutions and also the position where the smaller
  //map fits into to big one.
  //of course the cutposition should match the size of the fine map

  int coarsex = coarsemap->size2;
  int coarsey = coarsemap->size1;
  int finex = outmap->size2;
  int finey = outmap->size1;

  double xratio = (double) coarsex/finex;
  double yratio = (double) coarsey/finey;
  int coarseposx;
  int coarseposy;


  for(int i = 0; i < finey; i++)
    {
      for(int j = 0; j < finex; j++)
	{
	  coarseposy = (int) floor(i*yratio);
	  coarseposx = (int) floor(j*xratio);
	  if((i >= y1pos && i < y2pos) && (j >= x1pos && j < x2pos))
	    {
	      gsl_matrix_set(outmap,i,j,gsl_matrix_get(finemap,i-y1pos,j-x1pos));
	    }
	  else
	    {
	      gsl_matrix_set(outmap,i,j,gsl_matrix_get(coarsemap,coarseposy,coarseposx));
	    }
	}
    }
}

double cosmicweight(double lensredshift, double sourceredshift)
{

  double cosmicweight;
  cosModel cosmo1(0.3,0.7,-1.0);
  cosmicweight = (cosmo1.angularDistance(lensredshift,sourceredshift)*cosmo1.angularDistance(0,100000.0))/(cosmo1.angularDistance(0,sourceredshift)*cosmo1.angularDistance(lensredshift,100000.0));

  return cosmicweight;
}

void ReadMsystemInfo(string input, gsl_matrix *output)
{

  ifstream input1(input.c_str());
  string line;
  double value1;
  int syscounter = -1;
  int index = 0;

  while(input1)
    {
      getline(input1,line);
      
      if(line.size() > 1)
	{
	  if(line.at(0) == '#')
	    {
	      syscounter++;
	      index = 0;
	    }
	  else
	    {
	      istringstream read1(line);
	      //istringstream read2(line.substr(line.find_first_of(',')+1,line.size()));
	      read1 >> value1;
	      //read2 >> value2;
	      gsl_matrix_set(output,index,syscounter,value1);
	      index++;
	      //gsl_matrix_set(output,index,syscounter,value2);
	      //index++;
	    }
	}
    }
  input1.close();
}

void WriteMsystemInfo(gsl_matrix *input,string output)
{
  ofstream output1(output.c_str());

  for(int i = 0 ; i < input->size2; i++)
    {

      output1 <<"#" <<i+1 <<endl;

      for(int j = 0; j < input->size1; j++)
	{
	  output1 <<showpoint <<gsl_matrix_get(input,j,i) <<endl;
	}
    }

  output1.close();

}

void read_doubles(string input, int number, vector<double> &output)
{
    string::size_type index1;

    for(int i = 0; i < number; i++)
	{
	    istringstream helper;
	    string subline;
	    double value;
	    index1 = input.find_first_of("1234567890.-",0);
	    if(index1 != string::npos)
		{
		    input.erase(0,index1);
		    subline = input.substr(0,input.find_first_not_of("1234567890.-",0));
		    helper.str(subline);
		    if(subline.size() < 1)
			{
			    throw invalid_argument("Invalid double in string");
			}
		    input.erase(0,subline.size());
		    subline.clear();
		    helper >>value;
		    output.push_back(value);
		}
	    else
		{
		    throw invalid_argument("Not enough doubles in string");
		}
	}
}

void read_ints(string input, int number, vector<int> &output)
{
    string::size_type index1;

    for(int i = 0; i < number; i++)
	{
	    istringstream helper;
	    string subline;
	    int value;
	    index1 = input.find_first_of("1234567890-",0);
	    if(index1 != string::npos)
		{
		    input.erase(0,index1);
		    subline = input.substr(0,input.find_first_not_of("1234567890-",0));
		    helper.str(subline);
		    if(subline.size() < 1)
			{
			    throw invalid_argument("Invalid double in string");
			}
		    input.erase(0,subline.size());
		    subline.clear();
		    helper >>value;
		    output.push_back(value);
		}
	    else
		{
		    throw invalid_argument("Not enough doubles in string");
		}
	}
}

void read_doubles(string input, vector<double> &output)
{
    string::size_type index1 = 0;

    for(int i = 0; index1 != string::npos; i++)
	{
	    istringstream helper;
	    string subline;
	    double value;
	    index1 = input.find_first_of("1234567890.+-e",0);

	    input.erase(0,index1);
	    subline = input.substr(0,input.find_first_not_of("1234567890.+-e",0));
	    helper.str(subline);
	    if(subline.size() < 1)
		{
		    throw invalid_argument("Invalid double in string");
		}
	    input.erase(0,subline.size());
	    subline.clear();
	    helper >>value;
	    output.push_back(value);
	    index1 = input.find_first_of("1234567890.+-e",0);		
	}
}

void read_ints(string input, vector<int> &output)
{
    string::size_type index1 = 0;

    for(int i = 0; index1 != string::npos; i++)
	{
	    istringstream helper;
	    string subline;
	    int value;
	    index1 = input.find_first_of("1234567890-",0);

	    input.erase(0,index1);
	    subline = input.substr(0,input.find_first_not_of("1234567890-",0));
	    helper.str(subline);
	    if(subline.size() < 1)
		{
		    throw invalid_argument("Invalid integer in string");
		}
	    input.erase(0,subline.size());
	    subline.clear();
	    helper >>value;
	    output.push_back(value);
	    index1 = input.find_first_of("1234567890-",0);		
	}
}

double read_double(string input)
{
    string::size_type index1;
    istringstream helper;
    string subline;
    double value;
    index1 = input.find_first_of("1234567890.+-e",0);
    
    input.erase(0,index1);
    subline = input.substr(0,input.find_first_not_of("1234567890.+-e",0));
    helper.str(subline);
    if(subline.size() < 1)
	{
	    throw invalid_argument("Double in string invalid");
	}
    subline.clear();
    helper >>value;

    return value;
}

int read_int(string input)
{
    string::size_type index1;
    istringstream helper;
    string subline;
    int value;
    index1 = input.find_first_of("1234567890-",0);
    
    input.erase(0,index1);
    subline = input.substr(0,input.find_first_not_of("1234567890-",0));
    helper.str(subline);
    if(subline.size() < 1)
	{
	    throw invalid_argument("Integer in string invalid");
	}
    subline.clear();
    helper >>value;

    return value;
}



string read_word(string input, string symbol)
{

    string::size_type index1;
    string subline;

    index1 = input.find(symbol,0);
    if(index1 == string::npos)
	{
	    throw invalid_argument("Marker symbol not found in string");
	}
    input.erase(0,index1+symbol.size());
    subline = input.substr(0,input.find(symbol,0));
    
    return subline;
}

bool read_mind(string input)
{

    if((input.find("YES",0) != string::npos || input.find("yes",0) != string::npos || input.find("Jawohl",0) != string::npos || input.find("Yes",0) != string::npos || input.find("TRUE",0) != string::npos || input.find("true",0) != string::npos) && input.find("NO",0) == string::npos && input.find("no",0) == string::npos && input.find("FALSE",0) == string::npos && input.find("false",0) == string::npos && input.find("No",0) == string::npos && input.find("Nicht doch",0) == string::npos)
	{

	    return true;
	}
    else if(input.find_first_of("yYtT",0) < input.find_first_of("nNfF",0))
	{
	    return true;
	}
    else
	{
	    return false;
	}
}

void read_sequence(string input, vector<int> &output)
{

    string subline;
    output.clear();
    int value1, value2;

    while(input.size() != 0)
	{
	    subline = input.substr(0,input.find_first_of(","));
	    if(input.find(",",0) != string::npos)
		{
		    input.erase(0,input.find_first_of(",")+1);
		}
	    else
		{
		    input.clear();
		}


	    if(subline.find("-",0) != string::npos)
		{
		    istringstream helper1;
		    istringstream helper2;

		    helper1.str(subline.substr(0,subline.find_first_of("-",0)));
		    helper2.str(subline.substr(subline.find_first_of("-")+1));
		    helper1 >>value1;
		    helper2 >>value2;


		    for(int i = value1; i <= value2; i++)
			{
			    output.push_back(i);
			}  
		}
	    else
		{
		    istringstream helper1;
		    helper1.str(subline);
		    helper1 >>value1;
		    output.push_back(value1);
		}
	}
}



			       
   






