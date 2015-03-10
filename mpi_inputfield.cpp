#include "mpi_inputfield.h"

using namespace std;

string parametercut2(string input)
{
  int start;
  int stop;
  start = input.find("\"")+1;
  stop = input.rfind("\"")+1;
  string parameter;


  parameter.insert(0,input,start,stop-start-1);

  return parameter;
}

int countsymbol2(string input, const char &symbol)
{
  int number = 0;

  for(int i = 0; i < input.length(); i++)
    {
      if(input[i] == symbol)
	{
	  number++;
	}
    }
  return number;
}

void cutints2(string input,string symbol, gsl_vector_int *output, int number)
{
  string parameter = "";
  int cut;


  for(int i = 0; i < number-1; i++)
    {

      parameter.insert(0,input,0,input.find(symbol));
      input.erase(0,parameter.length()+1);
      istringstream para(parameter);
      para >>cut;
      gsl_vector_int_set(output,i,cut);
      parameter.clear();
      cut = 0;
    }
  istringstream para(input);
  para >>cut;
  gsl_vector_int_set(output,number-1,cut);
}

void cutdoubles2(string input,string symbol, gsl_matrix *output,int col, int number)
{
  string parameter = "";
  double cut;


  for(int i = 0; i < number-1; i++)
    {

      parameter.insert(0,input,0,input.find(symbol));
      input.erase(0,parameter.length()+1);
      istringstream para(parameter);
      para >>cut;
      gsl_matrix_set(output,i,col,cut);
      parameter.clear();
      cut = 0;
    }
  istringstream para(input);
  para >>cut;
  gsl_matrix_set(output,number-1,col,cut);
}


DataGrid::DataGrid(FieldOptions &options, int iterationindex,int my_rank,int p)
{

  //For destructor reasons:

  desshear = options.showbool("shear");
  desflexion = options.showbool("flexion");
  desccurve = options.showbool("ccurve");
  desmsystems = options.showbool("msystems");

  if(my_rank == 0)
    {
      root = true;
    }
  else
    {
      root = false;
    }


  //Calculating basic quantities

  fieldx = options.showdouble("fieldx2")-options.showdouble("fieldx1");
  x_dim = options.showstep(iterationindex);
  pixelsize = fieldx/x_dim;

  if(options.showbool("square"))
    {
      y_dim = x_dim;
    }
  else
    {
      y_dim = (int) ceil((options.showdouble("fieldy2")-options.showdouble("fieldy1"))/pixelsize);
    }
  fieldy = y_dim * pixelsize;
  fieldsize = fieldx * fieldy;
  xcentre = options.showdouble("fieldx1") + pixelsize/2.0;
  ycentre = options.showdouble("fieldy1") + pixelsize/2.0;
  maskcheck = gsl_matrix_int_calloc(y_dim,x_dim);



  // Checking the masking of the field, this is done right now because of memory allocation later.
  if(my_rank == 0)
    {
      string maskfile = options.showstring("maskinput");
      
      ifstream parafile(maskfile.c_str());
      string paraline;
      string parameter;
      bool badgrid;
      
      int rect;
      gsl_matrix *rectmasks;
      
      getline(parafile,paraline);
      istringstream par2(parametercut2(paraline));
      par2 >>rect;
      par2.str("");
      
      
      if(rect != 0)
	{
	  rectmasks = gsl_matrix_calloc(4,rect);
	  
	  for(int i = 0; i < rect; i++)
	    {
	      getline(parafile,paraline);
	      parameter = parametercut2(paraline);
	      cutdoubles2(parameter,",",rectmasks,i,4);
	    }
	}
      
      
      
      // Checking for each pixel centre if it is masked
      
      double xcoord;
      double ycoord;
      bool check;
      fieldpixels = 0;
      
      
      for(int i = 0; i < y_dim; i++)
	{
	  for(int j = 0; j < x_dim; j++)
	    {
	      check = false;
	      xcoord = (double) xcentre + (double) j*pixelsize;
	      ycoord = (double) ycentre + (double) i*pixelsize;
	      
	      
	      
	      if(rect != 0)
		{
		  for(int k = 0; k < rect; k++)
		    {
		      if((xcoord >= gsl_matrix_get(rectmasks,0,k)) && (xcoord <= gsl_matrix_get(rectmasks,1,k)) && (ycoord >= gsl_matrix_get(rectmasks,2,k)) && (ycoord <= gsl_matrix_get(rectmasks,3,k)))
			{
			  check= true;
			}
		    }
		}
	      
	      //write result
	      
	      
	      if(check)
		{
		  gsl_matrix_int_set(maskcheck,i,j,1);
		}
	      else
		{
		  gsl_matrix_int_set(maskcheck,i,j,0);
		}
	    }
	}
      
      if(rect != 0)
	{
	  gsl_matrix_free(rectmasks);
	} 
    

      //Trying to fit a findif grid
      FinDifGrid grid1(maskcheck,1.0,true);
      grid1.givegridmask(maskcheck);

      
      
      //converting the fieldmask
      
      for(int i = 0; i < y_dim; i++)
	{
	  for(int j = 0; j < x_dim; j++)
	    {
	      if(gsl_matrix_int_get(maskcheck,i,j) >= 0)
		{
		  gsl_matrix_int_set(maskcheck,i,j,1001+fieldpixels);
		  fieldpixels++;
		}
	      else
		{
		  gsl_matrix_int_set(maskcheck,i,j,1);
		}
	    }
	  
	}


  


      cout <<"Number of non-masked pixels: " <<fieldpixels <<endl;
    }

  //Sending the mask to all processes
  MPI_Barrier(MPI_COMM_WORLD);
  send_gsl_toworld(maskcheck);
  MPI_Barrier(MPI_COMM_WORLD);
  


  //reading the data and allocate memory for later analysis

  if(my_rank == 0)
    {
      if(options.showbool("shear"))
	{
	    string shearfile = options.show_filename("ellipinput",iterationindex);
	    ifstream shearinput(shearfile.c_str());
	  int shearcounter = 0;
	  double value1;
	  double value2;
	  double value3; 
	  double value4;
	  double value5;
	  while(shearinput >>value1 >>value2 >>value3 >>value4 >>value5)
	    {
	      if(value1 >= options.showdouble("areax1") && value1 <= options.showdouble("areax2") && value2 >= options.showdouble("areay1") && value2 <= options.showdouble("areay2"))
		{
		  shearcounter++;
		}
	    }
	  cout <<shearcounter <<"  ellipticity values read." <<endl;
	  numgalshear = shearcounter;
	  galaxydensityshear = (double) numgalshear/fieldsize;
	  shearinput.close();
	  rec = gsl_vector_calloc(numgalshear);
	  dec = gsl_vector_calloc(numgalshear);
	  ellip1 = gsl_vector_calloc(numgalshear);
	  ellip2 = gsl_vector_calloc(numgalshear);
	  shearweight = gsl_vector_calloc(numgalshear);
	  
	  ifstream shearread(shearfile.c_str());
	  value1 = 0.0;
	  value2 = 0.0;
	  value3 = 0.0; 
	  value4 = 0.0;
	  value5 = 0.0;
	  shearcounter = 0;
	  
	  while(shearread >>value1 >>value2 >>value3 >>value4 >>value5)
	    {
	      if(value1 >= options.showdouble("areax1") && value1 <= options.showdouble("areax2") && value2 >= options.showdouble("areay1") && value2 <= options.showdouble("areay2"))  
		{
		  gsl_vector_set(rec,shearcounter,value1);
		  gsl_vector_set(dec,shearcounter,value2);
		  gsl_vector_set(ellip1,shearcounter,value3);
		  gsl_vector_set(ellip2,shearcounter,value4);
		  gsl_vector_set(shearweight,shearcounter,value5);
		  shearcounter++;
		}
	    }
	  shearread.close();	  
	  
	}
      
      if(options.showbool("flexion"))
	{
	  
	  string flexionfile = options.show_filename("flexioninput",iterationindex);
	  ifstream flexioninput(flexionfile.c_str());
	  int flexioncounter = 0;
	  double value1;
	  double value2;
	  double value3; 
	  double value4;
	  double value5;
	  double value6;
	  double value7;
	  while(flexioninput >>value1 >>value2 >>value3 >>value4 >>value5 >>value6 >>value7)
	    {
	      if(value1 >= options.showdouble("areax1") && value1 <= options.showdouble("areax2") && value2 >= options.showdouble("areay1") && value2 <= options.showdouble("areay2"))
		{
		  flexioncounter++;
		}
	    }
	  cout <<flexioncounter <<"  flexion values read." <<endl;
	  numgalflexion = flexioncounter;
	  galaxydensityflexion = (double) numgalflexion/fieldsize;
	  flexioninput.close();
	  recflexion = gsl_vector_calloc(numgalflexion);
	  decflexion = gsl_vector_calloc(numgalflexion);
	  f1 = gsl_vector_calloc(numgalflexion);
	  f2 = gsl_vector_calloc(numgalflexion);
	  g1 = gsl_vector_calloc(numgalflexion);
	  g2 = gsl_vector_calloc(numgalflexion);
	  flexionweight = gsl_vector_calloc(numgalflexion);
	  
	  ifstream flexionread(flexionfile.c_str());
	  value1 = 0.0;
	  value2 = 0.0;
	  value3 = 0.0; 
	  value4 = 0.0;
	  value5 = 0.0;
	  value6 = 0.0;
	  value7 = 0.0;
	  flexioncounter = 0;
	  
	  while(flexionread >>value1 >>value2 >>value3 >>value4 >>value5 >>value6 >>value7)
	    {
	      if(value1 >= options.showdouble("areax1") && value1 <= options.showdouble("areax2") && value2 >= options.showdouble("areay1") && value2 <= options.showdouble("areay2"))  
		{
		  gsl_vector_set(recflexion,flexioncounter,value1);
		  gsl_vector_set(decflexion,flexioncounter,value2);
		  gsl_vector_set(f1,flexioncounter,value3);
		  gsl_vector_set(f2,flexioncounter,value4);
		  gsl_vector_set(g1,flexioncounter,value5);
		  gsl_vector_set(g2,flexioncounter,value6);
		  gsl_vector_set(flexionweight,flexioncounter,value7);
		  flexioncounter++;
		}
	    }
	  flexionread.close();	 	  
	}
      
      if(options.showbool("ccurve"))
	{
	  
	  string ccurvefilename = options.show_filename("ccurveinput",iterationindex);
	  ifstream ccurveinput(ccurvefilename.c_str());
	  double value1;
	  double value2;
	  double value3;
	  double value4;
	  double value5;
	  double value6;
	  int ccurvecounter = 0;
	  
	  //input reading
	  
	  while(ccurveinput >>value1 >>value2 >>value3 >>value4 >>value5 >>value6)
	    {
	      ccurvecounter++;
	    }
	  numccurvepts = ccurvecounter;
	  cout <<ccurvecounter <<"  Critical curve estimator points read." <<endl;
	  ccurveinput.close();
	  
	  ifstream ccurveread(ccurvefilename.c_str());
	  ccurverec = gsl_vector_calloc(numccurvepts);
	  ccurvedec = gsl_vector_calloc(numccurvepts);
	  ccurvered = gsl_vector_calloc(numccurvepts);
	  
	  value1 = 0.0;
	  value2 = 0.0;
	  value3 = 0.0;
	  value4 = 0.0;
	  value5 = 0.0;
	  value6 = 0.0;
	  
	  int index = 0;
	  while(ccurveread >>value1 >>value2 >>value3 >>value4 >>value5 >>value6)
	    {
	      gsl_vector_set(ccurverec,index,value1);
	      gsl_vector_set(ccurvedec,index,value2);
	      gsl_vector_set(ccurvered,index,value3);
	      index++;
	    }
	  
	  
	  ccurveread.close();
	  
    }
      
      if(options.showbool("msystems"))
	{
	  nummsystems = options.showint("nummsystems");
	  msysteminfo = gsl_matrix_calloc(18,options.showint("nummsystems"));
	  string msystemin = options.show_filename("msystemsinputinput",iterationindex);
	  int syscounter = -1;
	  int msindex = 0;
	  string line;
	  double value1;
	  double value2;
	  
	  ifstream msysteminput(msystemin.c_str());
	  
	  
	  while(msysteminput)
	    {
	      getline(msysteminput,line);
	      
	      if(line.size() > 1)
		{
		  if(line.at(0) == '#')
		    {
		      syscounter++;
		      msindex = 0;
		    }
		  else
		    {
		      istringstream read1(line.substr(0,line.find_first_of(',')));
		      istringstream read2(line.substr(line.find_first_of(',')+1,line.size()));
		      read1 >> value1;
		      read2 >> value2;
		      gsl_matrix_set(msysteminfo,msindex,syscounter,value1);
		      msindex++;
		      gsl_matrix_set(msysteminfo,msindex,syscounter,value2);
		      msindex++;
		    }
		}
	    }
	  msysteminput.close();
	}
    }
  MPI_Barrier(MPI_COMM_WORLD);


  //memory allocations and seind on all processes
  MPI_Bcast(&fieldpixels,1,MPI_INT,0,MPI_COMM_WORLD);


  if(options.showbool("shear"))
    {
      MPI_Bcast(&numgalshear,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Array_Manager smallarray(my_rank,p,y_dim,0);
      MPI_Array_Manager medarray(my_rank,p,fieldpixels,0);

      smallindex = smallarray.greenlight();
      medindex = medarray.greenlight();
      if(my_rank != 0)
	{
	  rec = gsl_vector_calloc(numgalshear);
	  dec = gsl_vector_calloc(numgalshear);
	  ellip1 = gsl_vector_calloc(numgalshear);
	  ellip2 = gsl_vector_calloc(numgalshear);
	  shearweight = gsl_vector_calloc(numgalshear);
	}
      if(smallarray.greenlight())
	{
	  meanellip1 = gsl_matrix_calloc(smallarray.showsize(),x_dim);
	  meanellip2 = gsl_matrix_calloc(smallarray.showsize(),x_dim); 
	  ellip1sd = gsl_matrix_calloc(smallarray.showsize(),x_dim);
	  ellip2sd = gsl_matrix_calloc(smallarray.showsize(),x_dim);
	  galaxyownersellip = gsl_matrix_uchar_calloc(smallarray.showsize()*x_dim, numgalshear);	  
	}
      else
	{
	  galaxyownersellip = gsl_matrix_uchar_calloc(1,1);
	}
      if(medarray.greenlight())
	{
	  galaxyshareellip = gsl_matrix_int_calloc(medarray.showsize(), fieldpixels);
	  ellip1covariance = gsl_matrix_calloc(medarray.showsize(), fieldpixels);
	  ellip2covariance = gsl_matrix_calloc(medarray.showsize(), fieldpixels);
	}
      internalellip1 = gsl_vector_calloc(fieldpixels);
      internalellip2 = gsl_vector_calloc(fieldpixels);
      internalellip1sd = gsl_vector_calloc(fieldpixels);
      internalellip2sd = gsl_vector_calloc(fieldpixels);
      finalgalaxyshareellip = gsl_matrix_int_calloc(fieldpixels,fieldpixels);
      finalgalaxyownersellip = gsl_matrix_uchar_calloc(x_dim*y_dim,numgalshear);



      if(my_rank == 0)
	{
	  finalmeanellip1 = gsl_matrix_calloc(y_dim,x_dim);
	  finalmeanellip2 = gsl_matrix_calloc(y_dim,x_dim);
	  finalellip1sd = gsl_matrix_calloc(y_dim,x_dim);
	  finalellip2sd = gsl_matrix_calloc(y_dim,x_dim);
	  finalellip1covariance = gsl_matrix_calloc(fieldpixels,fieldpixels);
	  finalellip2covariance = gsl_matrix_calloc(fieldpixels,fieldpixels);
	  pure_covariance1 = gsl_matrix_calloc(fieldpixels,fieldpixels);
	  pure_covariance2 = gsl_matrix_calloc(fieldpixels,fieldpixels);
	}
      else
	{
	  finalmeanellip1 = gsl_matrix_calloc(1,1);
	  finalmeanellip2 = gsl_matrix_calloc(1,1);
	  finalellip1sd = gsl_matrix_calloc(1,1);
	  finalellip2sd = gsl_matrix_calloc(1,1);
	  finalellip1covariance = gsl_matrix_calloc(1,1);
	  finalellip2covariance = gsl_matrix_calloc(1,1);
	  pure_covariance1 = gsl_matrix_calloc(1,1);
	  pure_covariance2 = gsl_matrix_calloc(1,1);
	}

	  
      
      send_gsl_toworld(rec);
      send_gsl_toworld(dec);
      send_gsl_toworld(ellip1);
      send_gsl_toworld(ellip2);
      send_gsl_toworld(shearweight);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  if(options.showbool("flexion"))
    {
      MPI_Bcast(&numgalflexion,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Array_Manager smallarray(my_rank,p,y_dim,0);
      MPI_Array_Manager medarray(my_rank,p,fieldpixels,0);

      smallindex = smallarray.greenlight();
      medindex = medarray.greenlight();

      if(my_rank != 0)
	{
	  recflexion = gsl_vector_calloc(numgalflexion);
	  decflexion = gsl_vector_calloc(numgalflexion);
	  f1 = gsl_vector_calloc(numgalflexion);
	  f2 = gsl_vector_calloc(numgalflexion);
	  g1 = gsl_vector_calloc(numgalflexion);
	  g2 = gsl_vector_calloc(numgalflexion);
	  flexionweight = gsl_vector_calloc(numgalflexion);
	}
      if(smallarray.greenlight())
	{
	  meanf1 = gsl_matrix_calloc(smallarray.showsize(), x_dim);
	  meanf2 = gsl_matrix_calloc(smallarray.showsize(), x_dim); 
	  f1sd = gsl_matrix_calloc(smallarray.showsize(), x_dim);
	  f2sd = gsl_matrix_calloc(smallarray.showsize(), x_dim);
	  meang1 = gsl_matrix_calloc(smallarray.showsize(), x_dim);
	  meang2 = gsl_matrix_calloc(smallarray.showsize(), x_dim); 
	  g1sd = gsl_matrix_calloc(smallarray.showsize(), x_dim);
	  g2sd = gsl_matrix_calloc(smallarray.showsize(), x_dim);
	  galaxyownersflexion = gsl_matrix_uchar_calloc(smallarray.showsize()*x_dim,numgalflexion);
	}
      else
	{
	  galaxyownersflexion = gsl_matrix_uchar_calloc(1,1);
	}


      if(medarray.greenlight())
	{
	  galaxyshareflexion = gsl_matrix_int_calloc(medarray.showsize(), fieldpixels);
	  f1covariance = gsl_matrix_calloc(medarray.showsize(), fieldpixels);
	  f2covariance = gsl_matrix_calloc(medarray.showsize(), fieldpixels);
	  g1covariance = gsl_matrix_calloc(medarray.showsize(), fieldpixels);
	  g2covariance = gsl_matrix_calloc(medarray.showsize(), fieldpixels);


	}
      internalf1 = gsl_vector_calloc(fieldpixels);
      internalf2 = gsl_vector_calloc(fieldpixels);
      internalf1sd = gsl_vector_calloc(fieldpixels);
      internalf2sd = gsl_vector_calloc(fieldpixels);
      internalg1 = gsl_vector_calloc(fieldpixels);
      internalg2 = gsl_vector_calloc(fieldpixels);
      internalg1sd = gsl_vector_calloc(fieldpixels);
      internalg2sd = gsl_vector_calloc(fieldpixels);
      finalgalaxyownersflexion = gsl_matrix_uchar_calloc(x_dim*y_dim,numgalflexion);
      finalgalaxyshareflexion = gsl_matrix_int_calloc(fieldpixels,fieldpixels);
      if(my_rank == 0)
	{
	  finalmeanf1 = gsl_matrix_calloc(y_dim,x_dim);
	  finalmeanf2 = gsl_matrix_calloc(y_dim,x_dim);
	  finalf1sd = gsl_matrix_calloc(y_dim,x_dim);
	  finalf2sd = gsl_matrix_calloc(y_dim,x_dim);
	  finalmeang1 = gsl_matrix_calloc(y_dim,x_dim);
	  finalmeang2 = gsl_matrix_calloc(y_dim,x_dim);
	  finalg1sd = gsl_matrix_calloc(y_dim,x_dim);
	  finalg2sd = gsl_matrix_calloc(y_dim,x_dim);
	  finalf1covariance = gsl_matrix_calloc(fieldpixels,fieldpixels);
	  finalf2covariance = gsl_matrix_calloc(fieldpixels,fieldpixels);
	  finalg1covariance = gsl_matrix_calloc(fieldpixels,fieldpixels);
	  finalg2covariance = gsl_matrix_calloc(fieldpixels,fieldpixels);
	}
      else
	{
	  finalmeanf1 = gsl_matrix_calloc(1,1);
	  finalmeanf2 = gsl_matrix_calloc(1,1);
	  finalf1sd = gsl_matrix_calloc(1,1);
	  finalf2sd = gsl_matrix_calloc(1,1);
	  finalmeang1 = gsl_matrix_calloc(1,1);
	  finalmeang2 = gsl_matrix_calloc(1,1);
	  finalg1sd = gsl_matrix_calloc(1,1);
	  finalg2sd = gsl_matrix_calloc(1,1);
	  finalf1covariance = gsl_matrix_calloc(1,1);
	  finalf2covariance = gsl_matrix_calloc(1,1);
	  finalg1covariance = gsl_matrix_calloc(1,1);
	  finalg2covariance = gsl_matrix_calloc(1,1);
	}

      
      send_gsl_toworld(recflexion);
      send_gsl_toworld(decflexion);
      send_gsl_toworld(f1);
      send_gsl_toworld(f2);
      send_gsl_toworld(g1);
      send_gsl_toworld(g2);
      send_gsl_toworld(flexionweight);
    }

  if(options.showbool("ccurve"))
    {
      MPI_Bcast(&numccurvepts,1,MPI_INT,0,MPI_COMM_WORLD);
      if(my_rank != 0)
	{
	  ccurverec = gsl_vector_calloc(numccurvepts);
	  ccurvedec = gsl_vector_calloc(numccurvepts);
	  ccurvered = gsl_vector_calloc(numccurvepts);
	}
      ccurve = gsl_matrix_int_calloc(y_dim, x_dim);
      ccurveerror = gsl_matrix_calloc(y_dim, x_dim);
      ccurveredshift = gsl_matrix_calloc(y_dim, x_dim);	
      send_gsl_toworld(ccurverec);
      send_gsl_toworld(ccurvedec);
      send_gsl_toworld(ccurvered);

    }
  if(options.showbool("msystems"))
    {
      MPI_Bcast(&nummsystems,1,MPI_INT,0,MPI_COMM_WORLD);
      if(my_rank != 0)
	{
	  msysteminfo = gsl_matrix_calloc(18,nummsystems);
	}
      send_gsl_toworld(msysteminfo);

    }



}

DataGrid::~DataGrid()
{

  //Deallocating memory, depending on the reconstruction type
  gsl_matrix_int_free(maskcheck);

  if(desshear)
    {
      gsl_vector_free(rec);
      gsl_vector_free(dec);
      gsl_vector_free(ellip1);
      gsl_vector_free(ellip2);
      gsl_vector_free(shearweight);
      gsl_matrix_uchar_free(galaxyownersellip);

      if(smallindex)
	{
	  gsl_matrix_free(meanellip1);
	  gsl_matrix_free(meanellip2);
	  gsl_matrix_free(ellip1sd);
	  gsl_matrix_free(ellip2sd);

	}
      gsl_vector_free(internalellip1);
      gsl_vector_free(internalellip2);
      gsl_vector_free(internalellip1sd);
      gsl_vector_free(internalellip2sd);
      gsl_matrix_uchar_free(finalgalaxyownersellip);
      gsl_matrix_int_free(finalgalaxyshareellip);
      if(medindex)
	{
	  gsl_matrix_int_free(galaxyshareellip);
	  gsl_matrix_free(ellip1covariance);
	  gsl_matrix_free(ellip2covariance);
	}
      gsl_matrix_free(finalmeanellip1);
      gsl_matrix_free(finalmeanellip2);
      gsl_matrix_free(finalellip1sd);
      gsl_matrix_free(finalellip2sd);
      gsl_matrix_free(finalellip1covariance);
      gsl_matrix_free(finalellip2covariance);
      gsl_matrix_free(pure_covariance1);
      gsl_matrix_free(pure_covariance2);
      
    }

  if(desflexion)
    {
      gsl_vector_free(recflexion);
      gsl_vector_free(decflexion);
      gsl_vector_free(f1);
      gsl_vector_free(f2);
      gsl_vector_free(g1);
      gsl_vector_free(g2);
      gsl_vector_free(flexionweight);
      if(smallindex)
	{
	  gsl_matrix_free(meanf1);
	  gsl_matrix_free(meanf2);
	  gsl_matrix_free(f1sd);
	  gsl_matrix_free(f2sd);
	  gsl_matrix_free(meang1);
	  gsl_matrix_free(meang2);
	  gsl_matrix_free(g1sd);
	  gsl_matrix_free(g2sd);      
	}
      gsl_vector_free(internalf1);
      gsl_vector_free(internalf2);
      gsl_vector_free(internalf1sd);
      gsl_vector_free(internalf2sd);
      gsl_vector_free(internalg1);
      gsl_vector_free(internalg2);
      gsl_vector_free(internalg1sd);
      gsl_vector_free(internalg2sd);
      gsl_matrix_uchar_free(galaxyownersflexion);
      if(medindex)
	{
	  gsl_matrix_int_free(galaxyshareflexion);
	  gsl_matrix_free(f1covariance);
	  gsl_matrix_free(f2covariance);
	  gsl_matrix_free(g1covariance);
	  gsl_matrix_free(g2covariance);
	}
      gsl_matrix_free(finalmeanf1);
      gsl_matrix_free(finalmeanf2);
      gsl_matrix_free(finalf1sd);
      gsl_matrix_free(finalf2sd);
      gsl_matrix_free(finalmeang1);
      gsl_matrix_free(finalmeang2);
      gsl_matrix_free(finalg1sd);
      gsl_matrix_free(finalg2sd);
      gsl_matrix_free(finalf1covariance);
      gsl_matrix_free(finalf2covariance);
      gsl_matrix_free(finalg1covariance);
      gsl_matrix_free(finalg2covariance);
      gsl_matrix_uchar_free(finalgalaxyownersflexion);
      gsl_matrix_int_free(finalgalaxyshareflexion);

      
    }


  if(desccurve)
    {
      gsl_vector_free(ccurverec);
      gsl_vector_free(ccurvedec); 
      gsl_vector_free(ccurvered);
      gsl_matrix_int_free(ccurve);
      gsl_matrix_free(ccurveerror);
      gsl_matrix_free(ccurveredshift);
    }

  if(desmsystems)
    {
      gsl_matrix_free(msysteminfo);
    }
    
}
void DataGrid::analyse(FieldOptions &options,int my_rank, int p)
{
  //int frac = (int) ceil((double) y_dim/p);
  if(options.showbool("shear"))
    {

      MPI_Array_Manager smallarray(my_rank,p,y_dim,0);
      MPI_Array_Manager medarray(my_rank,p,fieldpixels,0);
      smallarray.createcomm();
      medarray.createcomm();

      //Performing the averaging process first, afterwards the covariances are constructed

      //setting all output to the later masked value

      
      if(my_rank == 0)
	{
	  cout <<"Beginning adaptive shear averaging" <<endl;
	  cout <<"Progress:   |         | (of process 0)" <<endl;
	  CURSORUP(1);
	  CURSORRIGHT(13); 
	}
      //1st for loop over every pixel of your wanted grid
      

      double xcoord;
      double ycoord;
      double radius;
      int galcontained;

      int counter = 0;
      double progress = 0.1;

      for(int i = 0; i < smallarray.showsize(); i++)
	{
	  int ii = i + smallarray.showstartindex();
	  if(ii < y_dim)
	    {
	      if( (double) i/smallarray.showsize() > progress && my_rank == 0)
		{
		  cout <<"#" <<flush;
		  progress += 0.1;
		}
	      for(int j = 0; j < x_dim; j++)
		{
		  if(gsl_matrix_int_get(maskcheck,ii,j) != 1)
		    {
		      
		      //caluclate pixelcenter
		      
		      xcoord = xcentre + j*pixelsize;
		      ycoord = ycentre + ii*pixelsize;
		      galcontained = 0;
		      
		      //2nd for loop checking if enough galaxies contained
		      
		      for(int k = 0; galcontained < options.showint("sheargal"); k++)
			{
			  radius = (pixelsize/2.0)*(1.0+k*options.showdouble("radiusinc"));
			  galcontained = 0;
			  
			  //3rd loop over number of galaxies
			  for(int l = 0; l <numgalshear; l++)
			    {
			      if(sqrt(pow((xcoord-gsl_vector_get(rec,l)),2)+pow((ycoord-gsl_vector_get(dec,l)),2)) <= radius)
				{
				  galcontained++;
				}
			    }
			}
		      gsl_vector *ellip1sample = gsl_vector_calloc(galcontained);
		      gsl_vector *ellip2sample = gsl_vector_calloc(galcontained);
		      gsl_vector * weightsample =gsl_vector_calloc(galcontained);
		      int index = 0;
		      for(int k = 0; k < numgalshear; k++)
			{
			  
			  if(sqrt(pow((xcoord - gsl_vector_get(rec,k)),2)+pow((ycoord-gsl_vector_get(dec,k)),2)) <= radius)
			    {
			      gsl_vector_set(ellip1sample,index,gsl_vector_get(ellip1,k));
			      gsl_vector_set(ellip2sample,index,gsl_vector_get(ellip2,k));
			      gsl_vector_set(weightsample,index,gsl_vector_get(shearweight,k));
			      gsl_matrix_uchar_set(galaxyownersellip,counter,k,1);
			      index++;
			    }
			}
		      gsl_matrix_set(meanellip1,i,j,gsl_stats_wmean(gsl_vector_ptr(weightsample,0),1,gsl_vector_ptr(ellip1sample,0),1,index));
		      gsl_matrix_set(meanellip2,i,j,gsl_stats_wmean(gsl_vector_ptr(weightsample,0),1,gsl_vector_ptr(ellip2sample,0),1,index));
		      gsl_matrix_set(ellip1sd,i,j,gsl_stats_wsd(gsl_vector_ptr(weightsample,0),1,gsl_vector_ptr(ellip1sample,0),1,index));
		      gsl_matrix_set(ellip2sd,i,j,gsl_stats_wsd(gsl_vector_ptr(weightsample,0),1,gsl_vector_ptr(ellip2sample,0),1,index));		      
		    }
		  counter++;
		}
	    }
	}


      MPI_Barrier(MPI_COMM_WORLD);
      if(smallarray.greenlight())
	{
	  recv_gsl_fromworld(meanellip1,finalmeanellip1,smallarray);
	  recv_gsl_fromworld(meanellip2,finalmeanellip2,smallarray);
	  recv_gsl_fromworld(ellip1sd,finalellip1sd,smallarray);
	  recv_gsl_fromworld(ellip2sd,finalellip2sd,smallarray);
	}
      MPI_Barrier(MPI_COMM_WORLD);
      if(smallarray.greenlight())
	{
	  recv_gsl_fromworld(galaxyownersellip,finalgalaxyownersellip,smallarray);
	}
      send_gsl_toworld(finalgalaxyownersellip);
      MPI_Barrier(MPI_COMM_WORLD);



      if(my_rank == 0)
	{
	  int index = 0;
	  for(int i = 0; i < y_dim; i++)
	    {
	      for(int j = 0; j < x_dim; j++)
		{
		  if(gsl_matrix_int_get(maskcheck,i,j) != 1)
		    {
		      gsl_vector_set(internalellip1,index,gsl_matrix_get(finalmeanellip1,i,j));
		      gsl_vector_set(internalellip2,index,gsl_matrix_get(finalmeanellip2,i,j));
		      gsl_vector_set(internalellip1sd,index,gsl_matrix_get(finalellip1sd,i,j));
		      gsl_vector_set(internalellip2sd,index,gsl_matrix_get(finalellip2sd,i,j));
		      index++;
		    }
		}
	    }
	}

      int index = 0;
      int index2 = 0;

      for(int i = 0; i < y_dim; i++)
	{
	  for(int j = 0; j < x_dim; j++)
	    {
	      if(gsl_matrix_int_get(maskcheck,i,j) !=1)
		{
		  for(int k = 0; k < numgalshear; k++)
		    {
		      gsl_matrix_uchar_set(finalgalaxyownersellip,index,k,gsl_matrix_uchar_get(finalgalaxyownersellip,index2,k));
		    }
		  index++;
		}
	      index2++;
	    }
	}
		    
      MPI_Barrier(MPI_COMM_WORLD);
      send_gsl_toworld(internalellip1);
      send_gsl_toworld(internalellip2);
      send_gsl_toworld(internalellip1sd);
      send_gsl_toworld(internalellip2sd);
      MPI_Barrier(MPI_COMM_WORLD);


      if(my_rank == 0)
	{
	  cout <<endl;
	  cout <<"Performing ellipticity correlations: " <<endl;
	  cout <<"Progress:   |         | (of process 0)" <<endl;
	  CURSORUP(1);
	  CURSORRIGHT(13);
	} 

      int pos;
      progress = 0.1;
      for (int i = 0 ; i < medarray.showsize(); i++)
	{
	  int ii = i + medarray.showstartindex();
	  if(ii < fieldpixels)
	    {
	      if( (double) i/medarray.showsize() > progress && my_rank == 0)
		{
		  cout <<"#" <<flush;
		  progress += 0.1;
		}
	      
	      for(int j = -options.showint("pixelradius")*x_dim; j < options.showint("pixelradius")*x_dim; j++)
		{
		  
		  pos = ii + j;
		  int value = 0;
		  
		  if(pos < fieldpixels && pos >= ii)
		    {
		      for(int l = 0; l < numgalshear; l++)
			{
			  //if((gsl_matrix_int_get(galaxyownersellip,l,ii) == 1) && (gsl_matrix_int_get(galaxyownersellip,l,pos) == 1))
			  if((gsl_matrix_uchar_get(finalgalaxyownersellip,ii,l) == 1) && (gsl_matrix_uchar_get(finalgalaxyownersellip,pos,l) == 1))

			    {
			      value++;
			    }
			}		          
		      gsl_matrix_int_set(galaxyshareellip,i,pos,value);    

		    }
		}
	    }
	}
	  
      MPI_Barrier(MPI_COMM_WORLD);
      if(medarray.greenlight())
	{
	  recv_gsl_fromworld(galaxyshareellip,finalgalaxyshareellip,medarray);
	}
      if(my_rank == 0)
	{
	  for(int i =0; i < fieldpixels; i++)
	    {
	      for(int j = 0; j < i; j++)
		{
		  gsl_matrix_int_set(finalgalaxyshareellip,i,j,gsl_matrix_int_get(finalgalaxyshareellip,j,i));
		}
	    }
	}
      MPI_Barrier(MPI_COMM_WORLD);

      if(medarray.greenlight())
	{
	  send_gsl_toworld(finalgalaxyshareellip,medarray);
	}
      MPI_Barrier(MPI_COMM_WORLD);


	  for(int i = 0; i < medarray.showsize(); i++)  //calculating product of variance
	//weighted with number of shared galaxies
	{
	  int ii = i +  medarray.showstartindex();
	  if(i < fieldpixels)
	    {
	      for(int j = 0; j < fieldpixels; j++)
		{
		  gsl_matrix_set(ellip1covariance,i,j,(gsl_matrix_int_get(finalgalaxyshareellip,ii,j)/((gsl_matrix_int_get(finalgalaxyshareellip,ii,ii)+gsl_matrix_int_get(finalgalaxyshareellip,j,j))/2.0)*gsl_vector_get(internalellip1sd,ii)*gsl_vector_get(internalellip1sd,j)));
		  gsl_matrix_set(ellip2covariance,i,j,(gsl_matrix_int_get(finalgalaxyshareellip,ii,j)/((gsl_matrix_int_get(finalgalaxyshareellip,ii,ii)+gsl_matrix_int_get(finalgalaxyshareellip,j,j))/2.0)*gsl_vector_get(internalellip2sd,ii)*gsl_vector_get(internalellip2sd,j)));
		}
	    }
	
	}
      MPI_Barrier(MPI_COMM_WORLD);
	if(medarray.greenlight())
	  {
	    recv_gsl_fromworld(ellip1covariance,finalellip1covariance,medarray);
	    recv_gsl_fromworld(ellip2covariance,finalellip2covariance,medarray);
	  }
      MPI_Barrier(MPI_COMM_WORLD);
      if(my_rank == 0)
	{
	  cout <<endl;
	  cout <<"Inverting coefficient matrix..."<<flush;

/*
  ENTRY FOR CHARLES PHD PROJECT...NOT EFFICIENT!
*/

	  gsl_matrix_memcpy(pure_covariance1,finalellip1covariance);
	  gsl_matrix_memcpy(pure_covariance2,finalellip2covariance);
	  
	  invert_gsl(finalellip1covariance); //calling lapack routines in soph_math.h
	  invert_gsl(finalellip2covariance);
      
	  
	  cout <<"Done!" <<endl;
	}
      MPI_Barrier(MPI_COMM_WORLD);


    }


  if(options.showbool("flexion"))
    {
      //Performing the averaging process first, afterwards the covariances are constructed

      //setting all output to the later masked value

      MPI_Array_Manager smallarray(my_rank,p,y_dim,0);
      MPI_Array_Manager medarray(my_rank,p,fieldpixels,0);
      smallarray.createcomm();
      medarray.createcomm();

      if(my_rank == 0)
	{
 
	  cout <<"Beginning adaptive flexion averaging" <<endl;
	  cout <<"Progress:   |         | (of process 0)" <<endl;
	  CURSORUP(1);
	  CURSORRIGHT(13);
	} 

      //1st for loop over every pixel of your wanted grid
      double xcoord;
      double ycoord;
      double radius;
      int galcontained;
      int counter = 0;
      double progress = 0.1;

      for(int i = 0; i < smallarray.showsize(); i++)
	{
	  int ii = i + smallarray.showstartindex();
	  if(ii < y_dim)
	    {
	      if( (double) i/smallarray.showsize() > progress && my_rank == 0)
		{
		  cout <<"#" <<flush;
		  progress += 0.1;
		}
	      
	      for(int j = 0; j < x_dim; j++)
		{
		  if(gsl_matrix_int_get(maskcheck,ii,j) != 1)
		    {
		      
		      //caluclate pixelcenter
		      
		      xcoord = xcentre + j*pixelsize;
		      ycoord = ycentre + ii*pixelsize;
		      galcontained = 0;
		      
		      //2nd for loop checking if enough galaxies contained
		      
		      for(int k = 0; galcontained < options.showint("flexiongal"); k++)
			{
			  radius = (pixelsize/2.0)*(1.0+k*options.showdouble("radiusinc"));
			  galcontained = 0;
			  
			  //3rd loop over number of galaxies
			  for(int l = 0; l <numgalflexion; l++)
			    {
			      if(sqrt(pow((xcoord-gsl_vector_get(rec,l)),2)+pow((ycoord-gsl_vector_get(dec,l)),2)) <= radius)
				{
				  galcontained++;
				}
			    }
			}
		      gsl_vector *f1sample = gsl_vector_calloc(galcontained);
		      gsl_vector *f2sample = gsl_vector_calloc(galcontained);
		      gsl_vector *g1sample = gsl_vector_calloc(galcontained);
		      gsl_vector *g2sample = gsl_vector_calloc(galcontained);
		      gsl_vector * fweightsample =gsl_vector_calloc(galcontained);
		      int index = 0;
		      for(int k = 0; k < numgalflexion; k++)
			{
			  
			  if(sqrt(pow((xcoord - gsl_vector_get(rec,k)),2)+pow((ycoord-gsl_vector_get(dec,k)),2)) <= radius)
			    {
			      gsl_vector_set(f1sample,index,gsl_vector_get(f1,k));
			      gsl_vector_set(f2sample,index,gsl_vector_get(f2,k));
			      gsl_vector_set(g1sample,index,gsl_vector_get(g1,k));
			      gsl_vector_set(g2sample,index,gsl_vector_get(g2,k));
			      gsl_vector_set(fweightsample,index,gsl_vector_get(flexionweight,k));
			      gsl_matrix_uchar_set(galaxyownersflexion,counter,k,1);
			      index++;
			    }
			}
		      gsl_matrix_set(meanf1,i,j,gsl_stats_wmean(gsl_vector_ptr(fweightsample,0),1,gsl_vector_ptr(f1sample,0),1,index));
		      gsl_matrix_set(meanf2,i,j,gsl_stats_wmean(gsl_vector_ptr(fweightsample,0),1,gsl_vector_ptr(f2sample,0),1,index));
		      gsl_matrix_set(f1sd,i,j,gsl_stats_wsd(gsl_vector_ptr(fweightsample,0),1,gsl_vector_ptr(f1sample,0),1,index));
		      gsl_matrix_set(f2sd,i,j,gsl_stats_wsd(gsl_vector_ptr(fweightsample,0),1,gsl_vector_ptr(f2sample,0),1,index));
		      gsl_matrix_set(meang1,i,j,gsl_stats_wmean(gsl_vector_ptr(fweightsample,0),1,gsl_vector_ptr(g1sample,0),1,index));
		      gsl_matrix_set(meang2,i,j,gsl_stats_wmean(gsl_vector_ptr(fweightsample,0),1,gsl_vector_ptr(g2sample,0),1,index));
		      gsl_matrix_set(g1sd,i,j,gsl_stats_wsd(gsl_vector_ptr(fweightsample,0),1,gsl_vector_ptr(g1sample,0),1,index));
		      gsl_matrix_set(g2sd,i,j,gsl_stats_wsd(gsl_vector_ptr(fweightsample,0),1,gsl_vector_ptr(g2sample,0),1,index));
		    }
		  counter++;
		}
	    }
	}

      MPI_Barrier(MPI_COMM_WORLD);
      if(smallarray.greenlight())
	{
	  recv_gsl_fromworld(meanf1,finalmeanf1,smallarray);
	  recv_gsl_fromworld(meanf2,finalmeanf2,smallarray);
	  recv_gsl_fromworld(f1sd,finalf1sd,smallarray);
	  recv_gsl_fromworld(f2sd,finalf2sd,smallarray);
	  recv_gsl_fromworld(meang1,finalmeang1,smallarray);
	  recv_gsl_fromworld(meang2,finalmeang2,smallarray);
	  recv_gsl_fromworld(g1sd,finalg1sd,smallarray);
	  recv_gsl_fromworld(g2sd,finalg2sd,smallarray);
	}
      MPI_Barrier(MPI_COMM_WORLD);
      if(smallarray.greenlight())
	{
	  recv_gsl_fromworld(galaxyownersflexion,finalgalaxyownersflexion,smallarray);
	}
      send_gsl_toworld(finalgalaxyownersflexion);
      MPI_Barrier(MPI_COMM_WORLD);

      

      if(my_rank == 0)
	{
	  int index = 0;
	  for(int i = 0; i < y_dim; i++)
	    {
	      for(int j = 0; j < x_dim; j++)
		{
		  if(gsl_matrix_int_get(maskcheck,i,j) != 1)
		    {
		      gsl_vector_set(internalf1,index,gsl_matrix_get(finalmeanf1,i,j));
		      gsl_vector_set(internalf2,index,gsl_matrix_get(finalmeanf2,i,j));
		      gsl_vector_set(internalf1sd,index,gsl_matrix_get(finalf1sd,i,j));
		      gsl_vector_set(internalf2sd,index,gsl_matrix_get(finalf2sd,i,j));
		      gsl_vector_set(internalg1,index,gsl_matrix_get(finalmeang1,i,j));
		      gsl_vector_set(internalg2,index,gsl_matrix_get(finalmeang2,i,j));
		      gsl_vector_set(internalg1sd,index,gsl_matrix_get(finalg1sd,i,j));
		      gsl_vector_set(internalg2sd,index,gsl_matrix_get(finalg2sd,i,j));
		      index++;
		    }
		}
	    }
	}

      int index = 0;
      int index2 = 0;

      for(int i = 0; i < y_dim; i++)
	{
	  for(int j = 0; j < x_dim; j++)
	    {
	      if(gsl_matrix_int_get(maskcheck,i,j) !=1)
		{
		  for(int k = 0; k < numgalflexion; k++)
		    {
		      gsl_matrix_uchar_set(finalgalaxyownersflexion,index,k,gsl_matrix_uchar_get(finalgalaxyownersflexion,index2,k));
		    }
		  index++;
		}
	      index2++;
	    }
	}


      send_gsl_toworld(internalf1);
      send_gsl_toworld(internalf2);
      send_gsl_toworld(internalf1sd);
      send_gsl_toworld(internalf2sd);
      send_gsl_toworld(internalg1);
      send_gsl_toworld(internalg2);
      send_gsl_toworld(internalg1sd);
      send_gsl_toworld(internalg2sd);
      
      if(my_rank == 0)
	{
	  cout <<endl;
	  cout <<"Performing flexion correlations: " <<endl;
	  cout <<"Progress:   |         | (of process 0)" <<endl;
	  CURSORUP(1);
	  CURSORRIGHT(13);
	} 
      
      int pos;
      progress = 0.1;

      for (int i = 0; i < medarray.showsize(); i++)
	{
	  int ii = i + medarray.showstartindex(); 
	  if(ii < fieldpixels)
	    {
	      if( (double) i/medarray.showsize() > progress && my_rank == 0)
		{
		  cout <<"#" <<flush;
		  progress += 0.1;
		}
	      
	      for(int j = -options.showint("pixelradius")*x_dim; j < options.showint("pixelradius")*x_dim; j++)
		{
		  
		  pos = ii + j;
		  int value = 0;
		  
		  if(pos < fieldpixels && pos >= ii)
		    {
		      for(int l = 0; l < numgalflexion; l++)
			{
			  //if((gsl_matrix_int_get(galaxyownersflexion,l,ii) == 1) && (gsl_matrix_int_get(galaxyownersflexion,l,pos) == 1))
			  if((gsl_matrix_uchar_get(galaxyownersflexion,ii,l) == 1) && (gsl_matrix_uchar_get(galaxyownersflexion,pos,l) == 1))
			    {
			      value++;
			    }
			}
		      
		      gsl_matrix_int_set(galaxyshareflexion,i,pos,value);
		    }
		}
	    }
	}

      MPI_Barrier(MPI_COMM_WORLD);
      if(medarray.greenlight())
	{
	  recv_gsl_fromworld(galaxyshareflexion,finalgalaxyshareflexion,medarray);
	}
      for(int i = 0; i < fieldpixels; i++)
	{
	  for(int j = 0; j < i; j++)
	    {
	      gsl_matrix_int_set(finalgalaxyshareflexion,i,j,gsl_matrix_int_get(finalgalaxyshareflexion,j,i));
	    }
	   
	}
      MPI_Barrier(MPI_COMM_WORLD);
      if(medarray.greenlight())
	{
	  send_gsl_toworld(finalgalaxyshareflexion,medarray);
	}
      MPI_Barrier(MPI_COMM_WORLD);


      for(int i = 0; i < medarray.showsize(); i++)  //calculating product of variance
	//weighted with number of shared galaxies
	{
	  int ii = i + medarray.showstartindex();
	  if(ii < fieldpixels)
	    {
	      for(int j = 0; j < fieldpixels; j++)
		{
		  gsl_matrix_set(f1covariance,i,j,(gsl_matrix_int_get(finalgalaxyshareflexion,ii,j)/((gsl_matrix_int_get(finalgalaxyshareflexion,ii,ii)+gsl_matrix_int_get(finalgalaxyshareflexion,j,j))/2.0)*gsl_vector_get(internalf1sd,ii)*gsl_vector_get(internalf1sd,j)));
		  gsl_matrix_set(f2covariance,i,j,(gsl_matrix_int_get(finalgalaxyshareflexion,ii,j)/((gsl_matrix_int_get(finalgalaxyshareflexion,ii,ii)+gsl_matrix_int_get(finalgalaxyshareflexion,j,j))/2.0)*gsl_vector_get(internalf2sd,ii)*gsl_vector_get(internalf2sd,j)));
		  gsl_matrix_set(g1covariance,i,j,(gsl_matrix_int_get(finalgalaxyshareflexion,ii,j)/((gsl_matrix_int_get(finalgalaxyshareflexion,ii,ii)+gsl_matrix_int_get(finalgalaxyshareflexion,j,j))/2.0)*gsl_vector_get(internalg1sd,ii)*gsl_vector_get(internalg1sd,j)));
		  gsl_matrix_set(g2covariance,i,j,(gsl_matrix_int_get(finalgalaxyshareflexion,ii,j)/((gsl_matrix_int_get(finalgalaxyshareflexion,ii,ii)+gsl_matrix_int_get(finalgalaxyshareflexion,j,j))/2.0)*gsl_vector_get(internalg2sd,ii)*gsl_vector_get(internalg2sd,j)));
		}
	    }
	}

      MPI_Barrier(MPI_COMM_WORLD);
      if(medarray.greenlight())
	{
	  recv_gsl_fromworld(f1covariance,finalf1covariance,medarray);
	  recv_gsl_fromworld(f2covariance,finalf2covariance,medarray);
	  recv_gsl_fromworld(g1covariance,finalg1covariance,medarray);
	  recv_gsl_fromworld(g2covariance,finalg2covariance,medarray);
	}
      MPI_Barrier(MPI_COMM_WORLD);
      if(my_rank == 0)
	{
	  cout <<endl;
	  cout <<"Inverting coefficient matrix..."<<flush;
	  
	  invert_gsl(finalf1covariance); //calling lapack routines in soph_math.h
	  invert_gsl(finalf2covariance);
	  invert_gsl(finalg1covariance);
	  invert_gsl(finalg2covariance);
      
      
	  cout <<"Done!" <<endl;
	}


    }



  if(options.showbool("ccurve") && my_rank == 0)
    {

      gsl_matrix_int *gridpos = gsl_matrix_int_calloc(numccurvepts,2);
      int xpos = 0;
      int ypos = 0;
      int counter;


      //look for point position in field


      for(int i = 0; i < numccurvepts; i++)
	{
	  xpos = (int) floor((gsl_vector_get(ccurverec,i)-(xcentre-pixelsize/2.0))/pixelsize);
	  ypos = (int) floor((gsl_vector_get(ccurvedec,i)-(ycentre-pixelsize/2.0))/pixelsize);
	  gsl_matrix_int_set(gridpos,i,0,xpos);
	  gsl_matrix_int_set(gridpos,i,1,ypos);
	}


      for(int i = 0; i < y_dim; i++)
	{
	  for(int j = 0; j < x_dim; j++)
	    {
	      counter = 0;
	      for(int k = 0; k < numccurvepts; k++)
		{
		  if((gsl_matrix_int_get(gridpos,k,0)== j) && (gsl_matrix_int_get(gridpos,k,1)==i))
		    {
		      counter++;
		      gsl_matrix_set(ccurveredshift,i,j,gsl_matrix_get(ccurveredshift,i,j)+gsl_vector_get(ccurvered,k));
		    }
		}
	      if(gsl_matrix_int_get(maskcheck,i,j) != 1)
		{
		  gsl_matrix_int_set(ccurve,i,j,counter);
		  
		  if(counter != 0)
		    {
		      gsl_matrix_set(ccurveredshift,i,j,gsl_matrix_get(ccurveredshift,i,j)/counter);
		    }
		}
	    }
	}
      double cerror = pixelsize/30.0;  // see marcellos paper
      gsl_matrix_set_all(ccurveerror,cerror);


    }


  if(options.showbool("msystems") && my_rank == 0)
    {

      double xpos;
      double ypos;
      int counter;


      for(int i = 0; i < nummsystems; i++)
	{
	  gsl_matrix_set(msysteminfo,17,i,pixelsize);
	  for(int j = 0; j < (int) gsl_matrix_get(msysteminfo,0,i); j++)
	    {
	      counter = 0;
	      xpos = (double) floor((gsl_matrix_get(msysteminfo,2+2*j,i)-(xcentre-pixelsize/2.0))/pixelsize);
	      ypos = (double) floor((gsl_matrix_get(msysteminfo,3+2*j,i)-(ycentre-pixelsize/2.0))/pixelsize);
	      for(int k = 0; k < y_dim; k++)
		{
		  for(int l = 0; l < x_dim; l++)
		    {
		      if(gsl_matrix_int_get(maskcheck,k,l) != 1)
			{
			  counter++;
			}
		      if(k == ypos && l == xpos && gsl_matrix_int_get(maskcheck,k,l) != 1 )
			{
			  gsl_matrix_set(msysteminfo,12+j,i,(double) counter);
			}
		    }
		}
	    }
	}
      
    }
  MPI_Barrier(MPI_COMM_WORLD);
  
}

void DataGrid::write(FieldOptions &options, int iterationindex,int my_rank)
{
  if(my_rank == 0)
    {
      
      struct tm *date;
      time_t tim;
      time(&tim);
      
      date = localtime(&tim);
      string t = asctime(date);
      t = t.substr(0,24);
      
      
      if(options.showbool("shear"))
	{
	  
	  
	  cout <<"Writing ellipticity grids to FITS..."<<flush;
	  
//	  ostringstream output;
	  string outfile;
	  //  output <<options.showstring("shearoutput") <<options.showstep(iterationindex) <<".fits";
	  outfile = options.show_filename("shear",iterationindex);
	  
	  write_pimg(outfile,finalmeanellip1);
	  write_imge(outfile,"mean_ellip2",finalmeanellip2);
	  write_imge(outfile,"ellip1_sd",finalellip1sd);
	  write_imge(outfile,"ellip2_sd",finalellip2sd);
	  write_imgeint(outfile,"overlap",finalgalaxyshareellip);
	  write_imge(outfile,"ellip1_covariance",finalellip1covariance);
	  write_imge(outfile,"ellip2_covariance",finalellip2covariance);
	  write_imge(outfile,"ellip1_noninvcovariance",pure_covariance1);
	  write_imge(outfile,"ellip2_noninvcovariance",pure_covariance2);
	  write_imgeint(outfile,"field_mask",maskcheck);
	  
	  write_header(outfile,"DATE",t,"Creation time of file");
	  write_header(outfile,"OBJECT",options.showstring("clustername"),"Describing label of the reconstructed object");
	  write_header(outfile,"FILE",options.showstring("ellipinput"),"Ellipticity catalogue from which this file was created");
	  write_header(outfile,"MASK",options.showstring("maskinput"),"Mask definition used this file");
	  write_header(outfile,"FIELDX",fieldx,"Fieldsize in x-direction");
	  write_header(outfile,"FIELDY",fieldy,"Fieldsize in y-direction");
	  write_header(outfile,"AREA",fieldsize,"Area of the field");
	  write_header(outfile,"PIXELSIZE",pixelsize,"The sidelength of one pixel");
	  write_header(outfile,"PIXELS",fieldpixels,"Number of unmasked pixels in the field");
	  write_header(outfile,"EDGE_X",xcentre,"x coordinate of bottom-left pixel");
	  write_header(outfile,"EDGE_Y",ycentre,"y coordinate of bottom-left pixel");
	  write_header(outfile,"NUMGAL",numgalshear,"Number of galaxies in the field");
	  write_header(outfile,"GALDENS",galaxydensityshear,"Galaxy density in this field");
	  write_header(outfile,"GAL_PIXEL",options.showint("sheargal"),"Number of averaged galaxies per pixel");
	  
	  outfile = "";
	  
	  cout <<"Done." <<endl;
	}
      if(options.showbool("flexion"))
	{
	  
//	  ostringstream output;
	  string outfile;
	  //  output <<options.showstring("flexionoutput") <<options.showstep(iterationindex) <<".fits";
	  outfile = options.show_filename("flexion",iterationindex);
	  
	  cout <<"Writing flexion grids to FITS..." <<flush;
	  
	  write_pimg(outfile,finalmeanf1);
	  write_imge(outfile,"mean_f2",finalmeanf2);
	  write_imge(outfile,"f1_sd",finalf1sd);
	  write_imge(outfile,"f2_sd",finalf2sd);
	  write_imge(outfile,"mean_g1",finalmeang1);
	  write_imge(outfile,"mean_g2",finalmeang2);
	  write_imge(outfile,"g1_sd",finalg1sd);
	  write_imge(outfile,"g2_sd",finalg2sd);
	  write_imgeint(outfile,"overlap",finalgalaxyshareflexion);
	  write_imge(outfile,"f1_covariance",finalf1covariance);
	  write_imge(outfile,"f2_covariance",finalf2covariance);
	  write_imge(outfile,"g1_covariance",finalg1covariance);
	  write_imge(outfile,"g2_covariance",finalg2covariance);
	  write_imgeint(outfile,"field_mask",maskcheck);
	  write_header(outfile,"DATE",t,"Creation time of file");
	  write_header(outfile,"OBJECT",options.showstring("clustername"),"Describing label of the reconstructed object");
	  write_header(outfile,"FILE",options.showstring("flexioninput"),"Flexion catalogue from which this file was created");
	  write_header(outfile,"MASK",options.showstring("maskinput"),"Mask definition used this file");
	  write_header(outfile,"FIELDX",fieldx,"Fieldsize in x-direction");
	  write_header(outfile,"FIELDY",fieldy,"Fieldsize in y-direction");
	  write_header(outfile,"AREA",fieldsize,"Area of the field");
	  write_header(outfile,"PIXELSIZE",pixelsize,"The sidelength of one pixel");
	  write_header(outfile,"PIXELS",fieldpixels,"Number of unmasked pixels in the field");
	  write_header(outfile,"EDGE_X",xcentre,"x coordinate of bottom-left pixel");
	  write_header(outfile,"EDGE_Y",ycentre,"y coordinate of bottom-left pixel");
	  write_header(outfile,"NUMGAL",numgalflexion,"Number of galaxies in the field");
	  write_header(outfile,"GALDENS",galaxydensityflexion,"Galaxy density in this field");
	  write_header(outfile,"GAL_PIXEL",options.showint("flexiongal"),"Number of averaged galaxies per pixel");
	  
	  outfile = "";
	  
	  cout <<"Done." <<endl;
	}
      
      
      if(options.showbool("ccurve"))
	{
	  
//	  ostringstream output;
	  string outfile;
	  //  output <<options.showstring("ccurveoutput") <<options.showstep(iterationindex) <<".fits";
	  outfile = options.show_filename("ccurve",iterationindex);
	  
	  cout <<"Writing critical curve estimate grids to FITS..." <<flush;
	  
	  gsl_matrix_int *pccurve = gsl_matrix_int_calloc(y_dim,x_dim);
	  
	  for(int i = 0; i < y_dim; i++)
	    {
	      for(int j = 0; j < x_dim; j++)
		{
		  if(gsl_matrix_int_get(ccurve,i,j) != 0)
		    {
		      gsl_matrix_int_set(pccurve,i,j,1);
		    }
		  
		}
	    }
	  
	  write_pimgint(outfile,pccurve);
	  write_imgeint(outfile,"pixel_relevance",ccurve);
	  write_imge(outfile,"error",ccurveerror);
	  write_imge(outfile,"redshift_info",ccurveredshift);
	  write_imgeint(outfile,"field_mask",maskcheck);
	  write_header(outfile,"DATE",t,"Creation time of file");
	  write_header(outfile,"OBJECT",options.showstring("clustername"),"Describing label of the reconstructed object");
	  write_header(outfile,"FILE",options.showstring("ccurveinput"),"CCurve estimator file from which this file was created");
	  write_header(outfile,"MASK",options.showstring("maskinput"),"Mask definition used this file");
	  write_header(outfile,"FIELDX",fieldx,"Fieldsize in x-direction");
	  write_header(outfile,"FIELDY",fieldy,"Fieldsize in y-direction");
	  write_header(outfile,"AREA",fieldsize,"Area of the field");
	  write_header(outfile,"PIXELSIZE",pixelsize,"The sidelength of one pixel");
	  write_header(outfile,"PIXELS",fieldpixels,"Number of unmasked pixels in the field");
	  write_header(outfile,"EDGE_X",xcentre,"x coordinate of bottom-left pixel");
	  write_header(outfile,"EDGE_Y",ycentre,"y coordinate of bottom-left pixel");
	  write_header(outfile,"NUMCC",numccurvepts,"Number of points which are marked as part of critical curve");
	  write_header(outfile,"CCERROR",gsl_matrix_get(ccurveerror,0,0),"Error on the critical curve");
	  write_header(outfile,"MAXZ",gsl_matrix_max(ccurveredshift),"Highest redshift of a critical curve estimator");
	  
	  outfile = "";
	  cout <<"Done." <<endl;
	}
      
      if(options.showbool("msystems"))
	{
//	  ostringstream msoutput;
	  string fileout;
//	  msoutput <<options.showstring("msystemoutput") <<options.showstep(iterationindex) <<".dat";
	  fileout = options.show_filename("msystems",iterationindex);
	  cout <<"Writing multiple image systems info to ASCII..."<<flush;
	  
	  ofstream output(fileout.c_str());
	  
	  for(int i = 0 ; i < nummsystems; i++)
	    {
	      
	      output <<"#" <<i+1 <<endl;
	      
	      for(int j = 0; j < msysteminfo->size1; j++)
		{
		  output <<showpoint <<gsl_matrix_get(msysteminfo,j,i) <<endl;
		}
	      
	    }
	  
	  output.close();
	  fileout = "";
	  cout <<"Done." <<endl;
	  
	  
	}

    }
}














































