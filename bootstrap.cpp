#include "bootstrap.h"

Bootstrap::Bootstrap(int givenrealisationstart, int givenrealisationstop, gsl_matrix *input, string givenfilename)
{

  realisationstart = givenrealisationstart;
  realisationstop = givenrealisationstop;
  filename = givenfilename;
  length = input->size1;
  datalength = input->size2;
  data = gsl_matrix_calloc(length,datalength);
  gsl_matrix_memcpy(data,input);

}

Bootstrap::Bootstrap(int givenrealisationstart, int givenrealisationstop, CatWriter *input, string givenfilename)
{

  realisationstart = givenrealisationstart;
  realisationstop = givenrealisationstop;
  filename = givenfilename;
  input->create_copy();
  length = input->show_rowdim();
  datalength = input->show_coldim();
  data = gsl_matrix_calloc(length,datalength);
  //cout <<length <<"\t" <<datalength <<endl;
  //cout <<input->show_data()->size1 <<endl;
  gsl_matrix_memcpy(data,input->show_data());
  //cout <<"yay" <<endl;

}

Bootstrap::~Bootstrap()
{

  gsl_matrix_free(data);
  filename = "";

}


void Bootstrap::run(int sequencenumber)

{

  for(int i = realisationstart; i <= realisationstop; i++)
    {

      cout <<"Performing realisation: " <<i <<endl;

      int pos;
      time_t seed;
      time(&seed);
      const gsl_rng_type * T;
      gsl_rng * r;
      gsl_rng_env_setup();     
      T = gsl_rng_default;
      r = gsl_rng_alloc (T);
      gsl_rng_set(r,seed+i+sequencenumber); 
      

      ofstream output(filename.c_str());
      output <<scientific;
      output.width(12);


      for(int j = 0; j < length; j++)
	{

	  pos = (int) floor(gsl_rng_uniform(r)*length);
	  for(int l = 0; l < datalength; l++)
	      {

		  output <<gsl_matrix_get(data,pos,l) <<"\t" <<flush;
	      }
	  output <<endl;
	}

      output.close();

    }
}

void Bootstrap::run_delta(int sequencenumber)
{

    vector<int> BSsample;
    vector<int> goodsample;

    for(int i = 0; i < length; i++)
	{
	    if(gsl_matrix_get(data,i,5) == 1.0)
		{
		    BSsample.push_back(i);
		}
	    else
		{
		    goodsample.push_back(i);
		}
	}

    if((goodsample.size() + BSsample.size()) != length)
	{
	    throw runtime_error("BS samples are ill-defined.");
	}

    // cout <<goodsample.size() <<"\t" <<BSsample.size() <<endl;

    gsl_matrix *aux = gsl_matrix_calloc(length,6);
    
    for(int i = realisationstart; i <= realisationstop; i++)
	{
	    
	    cout <<"Performing realisation: " <<i <<endl;
	    
	    double pos;
	    int pos2;
	    time_t seed;
	    time(&seed);
	    const gsl_rng_type * T;
	    gsl_rng * r;
	    gsl_rng_env_setup();     
	    T = gsl_rng_default;
	    r = gsl_rng_alloc (T);
	    gsl_rng_set(r,seed+i+sequencenumber); 
	    
	    
	    ofstream output(filename.c_str());
	    output <<scientific;
	    output.width(12);
	    

	    for(int j = 0; j < length; j++)
		{
		    if(j == 0 || j == length -1 || gsl_matrix_get(data,j,4) != gsl_matrix_get(data,j-1,4))
			{
			    pos = (gsl_rng_uniform(r)*2.0 -1.0) * gsl_matrix_get(data,j,3);
			    pos += gsl_matrix_get(data,j,2);
			}
		    
		    for(int l = 0; l < 6; l++)
			{
			    if(l == 2)
				{
				    gsl_matrix_set(aux,j,l,pos); 
				}
			    else
				{
				    gsl_matrix_set(aux,j,l,gsl_matrix_get(data,j,l));
				}
			    
			}
		}
	    
	    for(int j = 0; j < length; j++)
		{
		    if(gsl_matrix_get(data,j,5) == 1.0)
			{
			    pos2 = BSsample[floor((double) gsl_rng_uniform(r)*BSsample.size())];
			}
		    else
			{
			    pos2 = j;
			}
		    for(int l = 0; l < 6; l++)
			{
			    output <<gsl_matrix_get(aux,pos2,l) <<"\t" <<flush;

			}
		    output <<endl;
		}
	    output.close();
	    
	}
}

void Bootstrap::set_length(int givenlength)
{
  length = givenlength;
}

void Bootstrap::set_datalength(int givendatalength)
{
  datalength = givendatalength;
}

void Bootstrap::set_realisations(int givenrealisationstart, int givenrealisationstop)
{
  realisationstart = givenrealisationstart;
  realisationstop = givenrealisationstop;
}

void Bootstrap::set_filename(string givenfilename)
{
  filename = givenfilename;
}

int Bootstrap::show_length()
{
  return length;
}

int Bootstrap::show_datalength()
{
  return datalength;
}

int Bootstrap::show_realisations(int selection)
{
  if(selection == 0)
    {
      return realisationstart;
    }
  else
    {
      return realisationstop;
    }

}

string Bootstrap::show_filename()
{
  return filename;
}
/*
BootstrapAnalysis::BootstrapAnalysis(int lengthgiven, int initresgiven, int interresgiven, string inputfilegiven)
{

  length = lengthgiven;
  initres = initresgiven;
  interres = interresgiven;
  inputfile = inputfilegiven;

  meanfield = gsl_vector_calloc(interres*interres);
  resultmap = gsl_vector_calloc(interres*interres);


  string filename;

  datamatrix = gsl_matrix_calloc(interres*interres,length);

  gsl_matrix *init = gsl_matrix_calloc(initres,initres);
  gsl_matrix *inter = gsl_matrix_calloc(interres,interres);
  gsl_vector *write = gsl_vector_calloc(interres*interres);

  cout <<"Reading the data..." <<endl;

  for(int i = 0; i < length; i++)
    {
      ostringstream file;
      file <<inputfile <<"/rec_realisation" <<i+1 <<"_" <<initres <<".fits";
      filename = file.str();
      read_imge(filename,"Convergence",init);

      if(initres == interres)
	{
	  gsl_matrix_memcpy(inter,init);
	}
      else
	{
	  cout <<"Interpolating..." <<flush;
	  cspline_mapinterpol_smooth(initres,initres,interres,interres,init,inter);
	  cout <<"Done!" <<endl;
	}


      matvec(inter,interres,interres,write);
      gsl_matrix_set_col(datamatrix,i,write);
      filename = "";
    }
  gsl_matrix_free(init);
  gsl_matrix_free(inter);
  gsl_vector_free(write);

}

BootstrapAnalysis::~BootstrapAnalysis()
{
  inputfile ="";
  gsl_matrix_free(datamatrix);
}


void BootstrapAnalysis::analyse(int outputselection)
{

  cout <<"Analysing the realisations..." <<flush;


  for(int i = 0; i < interres*interres; i++)
    {
      gsl_vector_set(meanfield,i,gsl_stats_mean(gsl_matrix_ptr(datamatrix,i,0),1,length));
          
      if(outputselection == 0)
	{
	  gsl_vector_set(resultmap,i,gsl_stats_mean(gsl_matrix_ptr(datamatrix,i,0),1,length));
	}
      else if(outputselection == 1)
	{
	  gsl_vector_set(resultmap,i,gsl_stats_sd(gsl_matrix_ptr(datamatrix,i,0),1,length));
	}
      else if(outputselection == 2)
	{
	  gsl_vector_set(resultmap,i,gsl_stats_variance(gsl_matrix_ptr(datamatrix,i,0),1,length));
	  } 
    }


  cout <<"Done!" <<endl;

}

void BootstrapAnalysis::write(gsl_vector *output)
{

  gsl_vector_memcpy(output,resultmap);
}


void BootstrapAnalysis::write(string output)
{

  write_pimg(output,interres,interres,resultmap);
  write_imge(output,"Mean",interres,interres,meanfield);

}

*/












 


 
