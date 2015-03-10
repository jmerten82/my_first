#include "options.h"

using namespace std;

string parametercut(string input)
{
  int start;
  int stop;
  start = input.find("\"")+1;
  stop = input.rfind("\"")+1;
  string parameter;


  parameter.insert(0,input,start,stop-start-1);

  return parameter;
}

int countsymbol(string input, const char &symbol)
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

void cutints(string input,string symbol, gsl_vector_int *output, int number)
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

void cutdoubles(string input,string symbol, gsl_vector *output, int number)
{
  string parameter = "";
  double cut;


  for(int i = 0; i < number-1; i++)
    {

      parameter.insert(0,input,0,input.find(symbol));
      input.erase(0,parameter.length()+1);
      istringstream para(parameter);
      para >>cut;
      gsl_vector_set(output,i,cut);
      parameter.clear();
      cut = 0;
    }
  istringstream para(input);
  para >>cut;
  gsl_vector_set(output,number-1,cut);
}




FieldOptions::FieldOptions(const string& filename)
{

  ifstream parafile (filename.c_str());
  string paraline;
  string parameter;


  getline(parafile,paraline);
  parameter = parametercut(paraline);
  square = parameter == "Yes";

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  shear = parameter == "Yes";

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  flexion = parameter == "Yes";

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  lowresccurve = parameter == "Yes";

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  msystems = parameter == "Yes";

  getline(parafile,paraline);
  istringstream ms(parametercut(paraline));
  ms >>nummsystems;
  ms.str("");

  getline(parafile,paraline);
  ellipinput = parametercut(paraline);

  getline(parafile,paraline);
  flexioninput = parametercut(paraline);

  getline(parafile,paraline);
  ccurveinput = parametercut(paraline);

  getline(parafile,paraline);
  msysteminput = parametercut(paraline);

  getline(parafile,paraline);
  maskinput = parametercut(paraline);

  getline(parafile,paraline);
  outputprefix = parametercut(paraline);

  getline(parafile,paraline);
  ccurveprefix = parametercut(paraline);

  getline(parafile,paraline);
  flexionprefix = parametercut(paraline);

  getline(parafile,paraline);
  msystemprefix = parametercut(paraline);

  getline(parafile,paraline);
  clustername = parametercut(paraline);


  getline(parafile,paraline);
  parameter = parametercut(paraline);
  cuts = gsl_vector_calloc(4);
  cutdoubles(parameter,",",cuts,4);

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  cuts2 = gsl_vector_calloc(4);
  cutdoubles(parameter,",",cuts2,4);

  getline(parafile,paraline);
  istringstream gal(parametercut(paraline));
  gal >>numgal;
  gal.str("");

  getline(parafile,paraline);
  istringstream flexiongal(parametercut(paraline));
  flexiongal >>flexionnumgal;
  flexiongal.str("");

  getline(parafile,paraline);
  istringstream start(parametercut(paraline));
  start >>startdim;
  start.str("");

  getline(parafile,paraline);
  istringstream stop(parametercut(paraline));
  stop >>stopdim;
  stop.str("");

  getline(parafile,paraline);
  istringstream step(parametercut(paraline));
  step >>stepsize;
  step.str("");

  getline(parafile,paraline);
  istringstream radius(parametercut(paraline));
  radius >>radiusinc;
  radius.str("");

  getline(parafile,paraline);
  istringstream covradius(parametercut(paraline));
  covradius >>covpixelradius;
  covradius.str("");

  getline(parafile,paraline);
  istringstream katze(parametercut(paraline));
  katze >>numbootstraps;
  katze.str("");

  parafile.close();

  steps = gsl_vector_int_calloc((int) ceil((double)(stopdim-startdim)/stepsize)+1);
  for(int i = 0; i < (int) ceil((double)(stopdim-startdim)/stepsize); i++)
    {
      gsl_vector_int_set(steps,i,startdim+i*stepsize);
    }
  gsl_vector_int_set(steps,(int) ceil((double)(stopdim-startdim)/stepsize),stopdim);

  bootstrapiterator = 0;
  


}

FieldOptions::~FieldOptions()
{

  ellipinput = "";
  flexioninput ="";
  ccurveinput = "";
  msysteminput = "";
  maskinput = "";
  outputprefix = "";
  flexionprefix ="";
  msystemprefix ="";
  gsl_vector_free(cuts);
  gsl_vector_free(cuts2);
  gsl_vector_int_free(steps);

}

void FieldOptions::setbool(string selection, bool value)
{
  if(selection == "shear")
    {
      shear = value;
    }
  else if(selection == "flexion")
    {
      flexion = value;
    }
  else if(selection == "msystems")
    {
      msystems = value;
    }
  else if(selection == "ccurve")
    {
      lowresccurve = value;
    }
  else if(selection == "square")
    {
      square = value;
    }
  else
    {
      throw invalid_argument("Selection not valid for setbool");
    }
}

bool FieldOptions::showbool(string selection)
{
  if(selection == "shear")
    {
      return shear;
    }
  else if(selection == "flexion")
    {
      return flexion;
    }
  else if(selection == "msystems")
    {
      return msystems;
    }
  else if(selection == "ccurve")
    {
      return lowresccurve;
    }
  else if(selection == "square")
    {
      return square;
    }
  else
    {
      throw invalid_argument("Selection not valid for showbool");
    }
}

void FieldOptions::setstring(string selection, string value)
{

  if(selection == "ellipinput")
    {
      ellipinput = value;
    }
  else if(selection == "flexioninput")
    {
      flexioninput = value;
    }
  else if(selection == "ccurveinput")
    {
      ccurveinput = value;
    }
  else if(selection == "msysteminput")
    {
      msysteminput = value;
    }
  else if(selection == "maskinput")
    {
      maskinput = value;
    }
  else if(selection == "shearoutput")
    {
      outputprefix = value;
    }
  else if(selection == "flexionoutput")
    {
      flexionprefix = value;
    }
  else if(selection == "ccurveoutput")
    {
      ccurveprefix = value;
    }
  else if(selection == "msystemoutput")
    {
      msystemprefix = value;
    }
  else if(selection == "clustername")
    {
      clustername = value;
    }
  else
    {
      throw invalid_argument("Selection not valid for setstring");
    }
}

string FieldOptions::showstring(string selection)
{
  if(selection == "ellipinput")
    {
      return ellipinput;
    }
  else if(selection == "flexioninput")
    {
      return flexioninput;
    }
  else if(selection == "ccurveinput")
    {
      return ccurveinput;
    }
  else if(selection == "msysteminput")
    {
      return msysteminput;
    }
  else if(selection == "maskinput")
    {
      return maskinput;
    }
  else if(selection == "shearoutput")
    {
      return outputprefix;
    }
  else if(selection == "flexionoutput")
    {
      return flexionprefix;
    }
  else if(selection == "ccurveoutput")
    {
      return ccurveprefix;
    }
  else if(selection == "msystemoutput")
    {
      return msystemprefix;
    }
  else if(selection == "clustername")
    {
      return clustername;
    }
  else
    {
      throw invalid_argument("Selection not valid for showstring");
    }
}

string FieldOptions::showstring(string selection,int iterationindex)
{
  if(selection == "ellipinput")
    {
      return ellipinput;
    }
  else if(selection == "flexioninput")
    {
      return flexioninput;
    }
  else if(selection == "ccurveinput")
    {
      return ccurveinput;
    }
  else if(selection == "msysteminput")
    {
      return msysteminput;
    }
  else if(selection == "maskinput")
    {
      return maskinput;
    }
  else if(selection == "shearoutput")
    {
      ostringstream output;
      output <<outputprefix <<gsl_vector_int_get(steps,iterationindex) <<".fits";
      return output.str();
    }
  else if(selection == "flexionoutput")
    {
      ostringstream output;
      output <<flexionprefix <<gsl_vector_int_get(steps,iterationindex) <<".fits";
      return output.str();
    }
  else if(selection == "ccurveoutput")
    {
      ostringstream output;
      output <<ccurveprefix <<gsl_vector_int_get(steps,iterationindex) <<".fits";
      return output.str();
    }
  else if(selection == "msystemoutput")
    {
      ostringstream output;
      output <<msystemprefix <<gsl_vector_int_get(steps,iterationindex) <<".fits";
      return output.str();
    }
  else if(selection == "clustername")
    {
      return clustername;
    }
  else
    {
      throw invalid_argument("Selection not valid for showstring");
    }
}

void FieldOptions::setdouble(string selection, double value)
{

  if(selection == "areax1")
    {
      gsl_vector_set(cuts,0,value);
    }
  else if(selection == "areax2")
    {
      gsl_vector_set(cuts,1,value);
    }
  else if(selection == "areay1")
    {
      gsl_vector_set(cuts,2,value);
    }
  else if(selection == "areay2")
    {
      gsl_vector_set(cuts,3,value);
    }
  else if(selection == "fieldx1")
    {
      gsl_vector_set(cuts2,0,value);
    }
  else if(selection == "fieldx2")
    {
      gsl_vector_set(cuts2,1,value);
    }
  else if(selection == "fieldy1")
    {
      gsl_vector_set(cuts2,2,value);
    }
  else if(selection == "fieldy2")
    {
      gsl_vector_set(cuts2,3,value);
    }
  else if(selection == "radiusinc")
    {
      radiusinc = value;
    }
  else
    {
      throw invalid_argument("Selection not valid for setdouble");
    }
}

double FieldOptions::showdouble(string selection)
{

  if(selection == "areax1")
    {
      return gsl_vector_get(cuts,0);
    }
  else if(selection == "areax2")
    {
      return gsl_vector_get(cuts,1);
    }
  else if(selection == "areay1")
    {
      return gsl_vector_get(cuts,2);
    }
  else if(selection == "areay2")
    {
      return gsl_vector_get(cuts,3);
    }
  else if(selection == "fieldx1")
    {
      return gsl_vector_get(cuts2,0);
    }
  else if(selection == "fieldx2")
    {
      return gsl_vector_get(cuts2,1);
    }
  else if(selection == "fieldy1")
    {
      return gsl_vector_get(cuts2,2);
    }
  else if(selection == "fieldy2")
    {
      return gsl_vector_get(cuts2,3);
    }
  else if(selection == "radiusinc")
    {
      return radiusinc;
    }
  else
    {
      throw invalid_argument("Selection not valid for showdouble");
    }
}

void FieldOptions::setint(string selection, int value)
{

  if(selection == "sheargal")
    {
      numgal = value;
    }
  else if(selection == "flexiongal")
    {
      flexionnumgal = value;
    }
  else if(selection == "startdim")
    {
      startdim = value;

      for(int i = 0; i < (int) ceil((double)(stopdim-startdim)/stepsize); i++)
	{
	  gsl_vector_int_set(steps,i,startdim+i*stepsize);
	}
      gsl_vector_int_set(steps,(int) ceil((double)(stopdim-startdim)/stepsize),stopdim);

    }
  else if(selection == "stopdim")
    {
      stopdim = value;

      for(int i = 0; i < (int) ceil((double)(stopdim-startdim)/stepsize); i++)
	{
	  gsl_vector_int_set(steps,i,startdim+i*stepsize);
	}
      gsl_vector_int_set(steps,(int) ceil((double)(stopdim-startdim)/stepsize),stopdim);

    }
  else if(selection == "stepsize")
    {
      stepsize = value;

      for(int i = 0; i < (int) ceil((double)(stopdim-startdim)/stepsize); i++)
	{
	  gsl_vector_int_set(steps,i,startdim+i*stepsize);
	}
      gsl_vector_int_set(steps,(int) ceil((double)(stopdim-startdim)/stepsize),stopdim);

    }
  else if(selection == "pixelradius")
    {
      covpixelradius = value;
    }
  else if( selection == "nummsystems")
    {
      nummsystems = value;
    }
  else if(selection == "bootstrap")
      {
	  numbootstraps = value;
      }
  else
    {
      throw invalid_argument("Selection not valid for setint");
    }
}

int FieldOptions::showint(string selection)
{
  if(selection == "sheargal")
    {
      return numgal;
    }
  else if(selection == "flexiongal")
    {
      return flexionnumgal;
    }
  else if(selection == "startdim")
    {
      return startdim;
    }
  else if(selection == "stopdim")
    {
      return stopdim;
    }
  else if(selection == "stepsize")
    {
      return stepsize;
    }
  else if(selection == "pixelradius")
    {
      return covpixelradius;
    }
  else if( selection == "nummsystems")
    {
      return nummsystems;
    }
  else if(selection == "bootstrap")
      {
	  return numbootstraps;
      }
  else if(selection == "BSiterator")
      {
	  return bootstrapiterator;
      }
  else
    {
      throw invalid_argument("Selection not valid for showint");
    }
}

int FieldOptions::showstep(int iterationindex)
{

  return gsl_vector_int_get(steps, iterationindex);
}

int FieldOptions::numres()
{

  return steps->size;
}

void FieldOptions::printall()
{

  cout <<"Here are all options on c_out: " <<endl;

  cout <<"Flags: " <<endl;
  cout <<"Shear is used: " <<shear <<endl;
  cout <<"Flexion is used: " <<flexion <<endl;
  cout <<"CCurve is used: " <<lowresccurve <<endl;
  cout <<"Multiple image systems are used: " <<msystems <<endl;
  cout <<"Field should be square: " <<square <<endl;

  cout <<endl;

  cout <<"Filenames: " <<endl;
  cout <<"Ellip cat: " <<ellipinput <<endl;
  cout <<"Flexion cat: " <<flexioninput <<endl;
  cout <<"CCurve input file: " <<ccurveinput <<endl;
  cout <<"MSytems input file: " <<msysteminput <<endl;
  cout <<"Mask input file: " <<maskinput <<endl;
  cout <<"Ellip field output: " <<outputprefix <<endl;
  cout <<"Flexion field output: " <<flexionprefix <<endl;
  cout <<"CCurve field output: " <<ccurveprefix <<endl;
  cout <<"MSytems output file: " <<msystemprefix <<endl;
  cout <<"Name of the object: " <<clustername <<endl;

  cout <<endl;

  cout <<"Double values: " <<endl;
  cout <<"Defining field points: " <<gsl_vector_get(cuts,0) <<" , " <<gsl_vector_get(cuts,1) <<" , " <<gsl_vector_get(cuts,2) <<" , " <<gsl_vector_get(cuts,3) <<endl;
  cout <<"Defining area points: " <<gsl_vector_get(cuts2,0) <<" , " <<gsl_vector_get(cuts2,1) <<" , " <<gsl_vector_get(cuts2,2) <<" , " <<gsl_vector_get(cuts2,3) <<endl;
  cout <<"Radius increment per step in adap. averaging: " <<radiusinc <<endl;

  cout <<endl;

  cout <<"Integer values: " <<endl;
  cout <<"Number of used shear galaxies per pixel: " <<numgal <<endl;
  cout <<"Number of used flexion galaxies per pixel: " <<flexionnumgal <<endl;
  cout <<"Start dimension: " <<startdim <<endl;
  cout <<"Stop dimension: " <<stopdim <<endl;
  cout <<"Stepsize between dimensions: " <<stepsize <<endl;
  cout <<"Radius to be searched while doing covariance: " <<covpixelradius <<endl;
  cout <<"Number of bootstrap realisations: "<<numbootstraps <<endl;

  cout <<endl;

  cout <<"Stepping: " <<endl;
  cout <<"Chosen resolutions: " <<flush;
  cout <<gsl_vector_int_get(steps,0) <<flush;
  for(int i = 1; i < steps->size; i++)
    {
      cout <<"," <<gsl_vector_int_get(steps,i) <<flush;
    }
  cout <<endl;
  cout <<"Number of chosen resolutions: " <<steps->size <<endl;
  cout <<"Number of multiple image systems: " <<nummsystems <<endl;

}

void FieldOptions::iterate()
{
    bootstrapiterator++;
}

string FieldOptions::show_filename(string selection, int olevel)
{

    string value;
    ostringstream helper;
    helper <<showstep(olevel);
    string helper2 = helper.str();
    ostringstream helper3;
    helper3 <<bootstrapiterator;
    string helper4 = helper3.str();

    if(selection == "shear")
	{
	    if(bootstrapiterator == 0)
		{
		    value = outputprefix+"_"+helper2+".fits";
		    return value;
		}
	    else
		{
		    value = outputprefix+"_"+helper2+"_BS"+helper4+".fits";
		    return value;
		}
	}

    else if(selection == "flexion")
	{
	    if(bootstrapiterator == 0)
		{
		    value = flexionprefix+"_"+helper2+".fits";
		    return value;
		}
	    else
		{
		    value = flexionprefix+"_"+helper2+"_BS"+helper4+".fits";
		    return value;
		}
	}
    else if(selection == "ccurve")
	{
	    if(bootstrapiterator == 0)
		{
		    value = ccurveprefix+"_"+helper2+".fits";
		    return value;
		}
	    else
		{
		    value = ccurveprefix+"_"+helper2+"_BS"+helper4+".fits";
		    return value;
		}
	}
    else if(selection == "msystems")
	{
	    if(bootstrapiterator == 0)
		{
		    value = msystemprefix+"_"+helper2+".fits";
		    return value;
		}
	    else
		{
		    value = msystemprefix+"_"+helper2+"_BS"+helper4+".fits";
		    return value;
		}
	}
    else if(selection == "ellipinput")
	{
	    if(bootstrapiterator == 0)
		{
		    value = ellipinput;
		    return value;
		}
	    else
		{
		    value = outputprefix+"_cat_BS"+helper4+".dat";
		    return value;
		}
	}
    else if(selection == "flexioninput")
	{
	    if(bootstrapiterator == 0)
		{
		    value = flexioninput;
		    return value;
		}
	    else
		{
		    value = flexionprefix+"_cat_BS"+helper4+".dat";
		    return value;
		}
	}
    else if(selection == "ccurveinput")
	{
	    if(bootstrapiterator == 0)
		{
		    value = ccurveinput;
		    return value;
		}
	    else
		{
		    value = ccurveprefix+"_cat_BS"+helper4+".dat";
		    return value;
		}
	}
    else if(selection == "msystemsinput")
	{
	    if(bootstrapiterator == 0)
		{
		    value = msysteminput;
		    return value;
		}
	    else
		{
		    value = msystemprefix+"_cat_BS"+helper4+".dat";
		    return value;
		}
	}
    else
	{
	    throw invalid_argument("Invalid selection for show_filename");
	}
}





ReconstructionOptions::ReconstructionOptions(const string& filename)
{
  ifstream parafile (filename.c_str());
  string paraline;
  string parameter;

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  shear = parameter == "Yes";

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  weaklimit = parameter == "Yes";
  
  getline(parafile,paraline);
  parameter = parametercut(paraline);
  flexion = parameter == "Yes";

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  ccurve = parameter == "Yes";

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  msystems = parameter == "Yes";

  getline(parafile,paraline);
  istringstream cut1(parametercut(paraline));
  cut1 >>nummsystems;
  cut1.str("");

  getline(parafile,paraline);
  istringstream cut2(parametercut(paraline));
  cut2 >>startdim;
  cut2.str("");

  getline(parafile,paraline);
  istringstream cut3(parametercut(paraline));
  cut3 >>stopdim;
  cut3.str("");

  getline(parafile,paraline);
  istringstream cut4(parametercut(paraline));
  cut4 >>stepsize;
  cut4.str("");

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  manualsteps = parameter == "Yes";

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  if(manualsteps)
    {
      steps = gsl_vector_int_calloc(countsymbol(parameter,',')+1);

      cutints(parameter,",",steps,countsymbol(parameter,',')+1);
      startdim = gsl_vector_int_get(steps,0);
      stopdim = gsl_vector_int_get(steps,steps->size -1);
    }
  else
    {
      steps = gsl_vector_int_calloc((int) ceil((double)(stopdim-startdim)/stepsize)+1);
      for(int i = 0; i < (int) ceil((double)(stopdim-startdim)/stepsize); i++)
	{
	  gsl_vector_int_set(steps,i,startdim+i*stepsize);
	}
      gsl_vector_int_set(steps,(int) ceil((double)(stopdim-startdim)/stepsize),stopdim);
    }

  getline(parafile,paraline);
  istringstream cut5(parametercut(paraline));
  cut5 >>iterations;
  cut5.str("");

  getline(parafile,paraline);
  istringstream cut6(parametercut(paraline));
  cut6 >>threshold;
  cut6.str("");

  getline(parafile,paraline);
  shearinput = parametercut(paraline);

  getline(parafile,paraline);
  flexioninput = parametercut(paraline);

  getline(parafile,paraline);
  ccurveinput = parametercut(paraline);

  getline(parafile,paraline);
  msysteminput = parametercut(paraline);

  getline(parafile,paraline);
  output = parametercut(paraline);

  getline(parafile,paraline);
  regtype = parametercut(paraline);

  getline(parafile,paraline);
  istringstream cut7(parametercut(paraline));
  cut7 >>regshear;
  cut7.str("");

  getline(parafile,paraline);
  istringstream cut8(parametercut(paraline));
  cut8 >>regflexion;
  cut8.str("");

  getline(parafile,paraline);
  istringstream cut9(parametercut(paraline));
  cut9 >>shearz;
  cut9.str("");

  getline(parafile,paraline);
  istringstream cut10(parametercut(paraline));
  cut10 >>flexionz;
  cut10.str("");

  getline(parafile,paraline);
  istringstream cut11(parametercut(paraline));
  cut11 >>clusterz;
  cut11.str("");

  getline(parafile,paraline);
  istringstream katze(parametercut(paraline));
  katze >>numbootstraps;
  katze.str("");

  bootstrapiterator = 0;
}

ReconstructionOptions::~ReconstructionOptions()
{
  gsl_vector_int_free(steps);
  shearinput ="";
  flexioninput ="";
  ccurveinput ="";
  msysteminput ="";
  output = "";
}

void ReconstructionOptions::setbool(string selection, bool value)
{
  if(selection == "shear")
    {
      shear = value;
    }

  else if(selection == "weaklimit")
    {
      weaklimit = value;
    }

  else if(selection == "flexion")
    {
      flexion = value;
    }
  else if(selection == "ccurve")
    {
      ccurve = value;
    }
  else if(selection == "msystems")
    {
      msystems = value;
    }
  else if(selection == "manualsteps")
    {
      manualsteps = value;
    }
  else
    {
      throw invalid_argument("Selection not valid for setbool");
    }
}

bool ReconstructionOptions::showbool(string selection)
{
  if(selection == "shear")
    {
      return shear;
    }

  else if(selection == "weaklimit")
    {
      return weaklimit;
    }

  else if(selection == "flexion")
    {
      return flexion;
    }
  else if(selection == "ccurve")
    {
      return ccurve;
    }
  else if(selection == "msystems")
    {
      return msystems;
    }
  else if(selection == "manualsteps")
    {
      return manualsteps;
    }
  else
    {
      throw invalid_argument("Selection not valid for showbool");
    }
}
 
void ReconstructionOptions::setdouble(string selection, double value)
{
  if(selection == "threshold")
    {
      threshold = value;
    }
  else if(selection == "regshear")
    {
      regshear = value;
    }
  else if(selection == "regflexion")
    {
      regflexion = value;
    }
  else if(selection == "shearz")
    {
      shearz = value;
    }
  else if(selection == "flexionz")
    {
      flexionz = value;
    }
  else if(selection == "clusterz")
    {
      clusterz = value;
    }
  else
    {
      throw invalid_argument("Selection not valid for setdouble");
    }
}

double ReconstructionOptions::showdouble(string selection)
{
  if(selection == "threshold")
    {
      return threshold;
    }
  else if(selection == "regshear")
    {
      return regshear;
    }
  else if(selection == "regflexion")
    {
      return regflexion;
    }
  else if(selection == "shearz")
    {
      return shearz;
    }
  else if(selection == "flexionz")
    {
      return flexionz;
    }
  else if(selection == "clusterz")
    {
      return clusterz;
    }
  else
    {
      throw invalid_argument("Selection not valid for showdouble");
    }
} 

void ReconstructionOptions::setint(string selection, int value)
{

  if(selection == "nummsystems")
    {
      nummsystems = value;
    }
  else if(selection == "startdim")
    {
      startdim = value;
    }
  else if(selection == "stopdim")
    {
      stopdim = value;
    }
  else if(selection == "stepsize")
    {
      stepsize = value;
    }
  else if(selection == "iterations")
    {
      iterations = value;
    }
  else if(selection == "bootstrap")
      {
	  numbootstraps = value;
      }
  else
    {
      throw invalid_argument("Selection not valid for setint");
    }
}

int ReconstructionOptions::showint(string selection)
{
  if(selection == "nummsystems")
    {
      return nummsystems;
    }
  else if(selection == "startdim")
    {
      return startdim;
    }
  else if(selection == "stopdim")
    {
      return stopdim;
    }
  else if(selection == "stepsize")
    {
      return stepsize;
    }
  else if(selection == "iterations")
    {
      return iterations;
    }
  else if(selection == "bootstrap")
      {
	  return numbootstraps;
      }
  else if(selection == "BSiterator")
      {
	  return bootstrapiterator;
      }
  else
    {
      throw invalid_argument("Selection not valid for showint");
    }
}

void ReconstructionOptions::setstring(string selection, string value)
{
  if(selection == "shearinput")
    {
      shearinput = value;
    }
  else if(selection == "flexioninput")
    {
      flexioninput = value;
    }
  else if(selection == "ccurveinput")
    {
      ccurveinput = value;
    }
  else if(selection == "msysteminput")
    {
      msysteminput = value;
    }
  else if(selection == "output")
    {
      output = value;
    }
  else if(selection == "regtype")
    {
      regtype = value;
    }
  else
    {
      throw invalid_argument("Selection not valid for setstring");
    }

}

string ReconstructionOptions::showstring(string selection)
{
  if(selection == "shearinput")
    {
      return shearinput;
    }
  else if(selection == "flexioninput")
    {
      return flexioninput;
    }
  else if(selection == "ccurveinput")
    {
      return ccurveinput;
    }
  else if(selection == "msysteminput")
    {
      return msysteminput;
    }
  else if(selection == "output")
    {
      return output;
    }
  else if(selection == "regtype")
    {
      return regtype;
    }
  else
    {
      throw invalid_argument("Selection not valid for showstring");
    }
}

string ReconstructionOptions::showstring(string selection, int iterationindex)
{
  if(selection == "shearinput")
    {
      ostringstream output;
      output <<shearinput <<gsl_vector_int_get(steps,iterationindex) <<".fits";
      return output.str();
    }
  else if(selection == "flexioninput")
    {
      ostringstream output;
      output <<flexioninput <<gsl_vector_int_get(steps,iterationindex) <<".fits";
      return output.str();
    }
  else if(selection == "ccurveinput")
    {
      ostringstream output;
      output <<ccurveinput <<gsl_vector_int_get(steps,iterationindex) <<".fits";
      return output.str();
    }
  else if(selection == "msysteminput")
    {
      ostringstream output;
      output <<msysteminput <<gsl_vector_int_get(steps,iterationindex) <<".fits";
      return output.str();
    }
  else if(selection == "output")
    {
      ostringstream output1;
      output1 <<output <<gsl_vector_int_get(steps,iterationindex) <<".fits";
      return output1.str();
    }
  else
    {
      throw invalid_argument("Selection not valid for showstring");
    }
}

void ReconstructionOptions::printall()
{

  cout <<"Here are all options on c_out: " <<endl;

  cout <<"Flags: " <<endl;
  cout <<"Shear is used: " <<shear <<endl;
  cout <<"Weak limit is used: " <<weaklimit <<endl;
  cout <<"Flexion is used: " <<flexion <<endl;
  cout <<"CCurve is used: " <<ccurve <<endl;
  cout <<"Multiple image systems are used: " <<msystems <<endl;
  cout <<"Steps are set manually: " <<manualsteps <<endl;

  cout <<endl;

  cout <<"Filenames: " <<endl;
  cout <<"Ellip cat: " <<shearinput <<endl;
  cout <<"Flexion cat: " <<flexioninput <<endl;
  cout <<"CCurve input file: " <<ccurveinput <<endl;
  cout <<"MSytems input file: " <<msysteminput <<endl;
  cout <<"Output files: " <<output <<endl;

  cout <<endl;

  cout <<"Double values: " <<endl;
  cout <<"Convergence change threshold: " <<threshold <<endl;
  cout <<"Shear regularisation parameter: " <<regshear <<endl;
  cout <<"Flexion regularisation parameter: " <<regflexion <<endl;
  cout <<"Average redshift of shear sources: " <<shearz <<endl;
  cout <<"Average redshift of flexion sources: " <<flexionz <<endl;
  cout <<"Redshift of the lens: " <<clusterz <<endl;

  cout <<endl;

  cout <<"Integer values: " <<endl;
  cout <<"Number of multiple image systems " <<nummsystems <<endl;
  cout <<"Start dimension: " <<startdim <<endl;
  cout <<"Stop dimension: " <<stopdim <<endl;
  cout <<"Maximum number of inner level iterations: " <<iterations <<endl;
  cout <<" Reconstruction steps: " <<flush;
  for(int i = 0; i < steps->size -1; i++)
    {
      cout <<gsl_vector_int_get(steps,i) <<"," <<flush;
    }
  cout <<gsl_vector_int_get(steps,steps->size -1) <<endl;

  cout <<"Number of bootstrap realisations: " <<endl;

}

int ReconstructionOptions::showstep(int iterationindex)
{
  return gsl_vector_int_get(steps,iterationindex);
}

int ReconstructionOptions::shownumsteps()
{
  return steps->size;
}

void ReconstructionOptions::iterate()
{

    bootstrapiterator++;
}

string ReconstructionOptions::show_filename(string selection, int olevel)
{

    string value;
    ostringstream helper;
    helper <<showstep(olevel);
    string helper2 = helper.str();
    ostringstream helper3;
    helper3 <<bootstrapiterator;
    string helper4 = helper3.str();

    if(selection == "shearinput")
	{
	    if(bootstrapiterator == 0)
		{
		    value = shearinput+"_"+helper2+".fits";
		    return value;
		}
	    else
		{
		    value = shearinput+"_"+helper2+"_BS"+helper4+".fits";
		    return value;
		}
	}
    else if(selection == "flexioninput")
	{
	    if(bootstrapiterator == 0)
		{
		    value = flexioninput+"_"+helper2+".fits";
		    return value;
		}
	    else
		{
		    value = flexioninput+"_"+helper2+"_BS"+helper4+".fits";
		    return value;
		}
	}

    else if(selection == "ccurveinput")
	{
	    if(bootstrapiterator == 0)
		{
		    value = ccurveinput+"_"+helper2+".fits";
		    return value;
		}
	    else
		{
		    value = ccurveinput+"_"+helper2+"_BS"+helper4+".fits";
		    return value;
		}
	}

    else if(selection == "msytemsinput")
	{
	    if(bootstrapiterator == 0)
		{
		    value = msysteminput+"_"+helper2+".fits";
		    return value;
		}
	    else
		{
		    value = msysteminput+"_"+helper2+"_BS"+helper4+".fits";
		    return value;
		}
	}
    else if(selection == "output")
	{
	    if(bootstrapiterator == 0)
		{
		    value = output+"_"+helper2+".fits";
		    return value;
		}
	    else
		{
		    value = output+"_"+helper2+"_BS"+helper4+".fits";
		    return value;
		}
	}

    else
	{
	    throw invalid_argument("Invalid selection for show_filename");
	}
}




HighresOptions::HighresOptions(const string& filename)
{
  ifstream parafile (filename.c_str());
  string paraline;
  string parameter;

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  ccurve = parameter == "Yes";
  
  getline(parafile,paraline);
  parameter = parametercut(paraline);
  msystems = parameter == "Yes";

  getline(parafile,paraline);
  lowresresult = parametercut(paraline);

  getline(parafile,paraline);
  highresccurve = parametercut(paraline);

  getline(parafile,paraline);
  highresmsystemsfile = parametercut(paraline);

  getline(parafile,paraline);
  interpolationoutput = parametercut(paraline);

  getline(parafile,paraline);
  msystemsoutput = parametercut(paraline);

  getline(parafile,paraline);
  mergedoutput = parametercut(paraline);

  getline(parafile,paraline);
  istringstream cut3(parametercut(paraline));
  cut3 >>finalx;
  cut3.str("");

  getline(parafile,paraline);
  istringstream cut4(parametercut(paraline));
  cut4 >>finaly;
  cut4.str("");

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  cuts = gsl_vector_int_calloc(4);
  cutints(parameter,",",cuts,4);

  getline(parafile,paraline);
  istringstream cut15(parametercut(paraline));
  cut15 >>nummsystems;
  cut15.str("");

  getline(parafile,paraline);
  istringstream cut5(parametercut(paraline));
  cut5 >>reg;
  cut5.str("");

  getline(parafile,paraline);
  regscheme = parametercut(paraline);

  getline(parafile,paraline);
  istringstream cut6(parametercut(paraline));
  cut6 >>clredshift;
  cut6.str("");

  getline(parafile,paraline);
  istringstream cut7(parametercut(paraline));
  cut7 >>ssigma;
  cut7.str("");

  getline(parafile,paraline);
  istringstream katze(parametercut(paraline));
  katze >>numbootstraps;
  katze.str("");

  bootstrapiterator = 0;

}

HighresOptions::~HighresOptions()
{
  gsl_vector_int_free(cuts);
  lowresresult = "";
  highresccurve = "";
  highresmsystemsfile = "";
  interpolationoutput = "";
  mergedoutput = "";
}

void HighresOptions::setbool(bool value,string selection)
{
  if(selection == "ccurve")
    {
      ccurve = value;
    }
  else if(selection == "msystems")
    {
      msystems = value;
    }
  else
    {
      throw invalid_argument("Selection not valid for setbool");
    }
}

void HighresOptions::setstring(string value,string selection)
{

  if(selection == "lowresresult")
    {
      lowresresult = value;
    }
  else if (selection == "highresccurve")
    {
      highresccurve = value;
    }
  else if(selection == "interpolationoutput")
    {
      interpolationoutput = value;
    }
  else if(selection == "mergedoutput")
    {
      mergedoutput = value;
    }
  else if (selection == "highresmsystemsfile")
    {
      highresmsystemsfile = value;
    }
  else if (selection == "msystemsoutput")
    {
      msystemsoutput = value;
    }
  else if (selection == "regscheme")
    {
      regscheme = value;
    }
  else
    {
      throw invalid_argument("Selection not valid for setstring");
    }

}

void HighresOptions::setint(int value,string selection)
{
  if(selection == "finalx")
    {
      finalx = value;
    }
  else if(selection == "finaly")
    {
      finaly = value;
    }
  else if(selection == "x_min")
    {
      gsl_vector_int_set(cuts,0,value);
    }
  else if(selection == "x_max")
    {
      gsl_vector_int_set(cuts,1,value);
    }
  else if(selection == "y_min")
    {
      gsl_vector_int_set(cuts,2,value);
    }
  else if(selection == "y_max")
    {
      gsl_vector_int_set(cuts,3,value);
    }
  else if(selection == "nummsystems")
    {
      nummsystems = value;
    }
  else if(selection == "bootstrap")
      {
	  numbootstraps = value;
      }
  else
    {
      throw invalid_argument("Selection not valid for setint");
    }
}
     

void HighresOptions::setdouble(double value, string selection)
{
  if(selection == "reg")
    {
      reg = value;
    }
  else if(selection == "clredshift")
    {
      clredshift = value;
    }
  else if(selection == "ssigma")
    {
      ssigma = value;
    }
  else
    {
      throw invalid_argument("Selection not valid for setdouble");
    }

}

bool HighresOptions::showbool(string selection)
{
  if(selection == "ccurve")
    {
      return ccurve;
    }
  else if(selection == "msystems")
    {
      return msystems;
    }
  else
    {
      throw invalid_argument("Selection not valid for showbool");
    }
}

string HighresOptions::showstring(string selection)
{

  if(selection == "lowresresult")
    {
      return lowresresult;
    }
  else if (selection == "highresccurve")
    {
      return highresccurve;
    }
  else if(selection == "interpolationoutput")
    {
      return interpolationoutput;
    }
  else if(selection == "mergedoutput")
    {
      return mergedoutput;
    }
  else if (selection == "highresmsystemsfile")
    {
      return highresmsystemsfile;
    }
  else if (selection == "msystemsoutput")
    {
      return msystemsoutput;
    }
  else if (selection == "regscheme")
    {
      return regscheme;
    }
  else
    {
      throw invalid_argument("Selection not valid for showstring");
    }

}

int HighresOptions::showint(string selection)
{
  if(selection == "finalx")
    {
      return finalx;
    }
  else if(selection == "finaly")
    {
      return finaly;
    }
  else if(selection == "x_min")
    {
      return gsl_vector_int_get(cuts,0);
    }
  else if(selection == "x_max")
    {
      return gsl_vector_int_get(cuts,1);
    }
  else if(selection == "y_min")
    {
      return gsl_vector_int_get(cuts,2);
    }
  else if(selection == "y_max")
    {
      return gsl_vector_int_get(cuts,3);
    }
  else if(selection == "nummsystems")
    {
      return nummsystems;
    }
  else if(selection == "bootstrap")
      {
	  return numbootstraps;
      }
  else if(selection == "BSiterator")
      {
	  return bootstrapiterator;
      }
  else
    {
      throw invalid_argument("Selection not valid for showint");
    }
}
     

double HighresOptions::showdouble(string selection)
{
  if(selection == "reg")
    {
      return reg;
    }
  else if(selection == "clredshift")
    {
      return clredshift;
    }
  else if(selection == "ssigma")
    {
      return ssigma;
    }
  else
    {
      throw invalid_argument("Selection not valid for showdouble");
    }

}

string HighresOptions::show_filename(string selection)
{

    string value;
    ostringstream helper3;
    helper3 <<bootstrapiterator;
    string helper4 = helper3.str();

    if(selection == "lowres_in")
	{

	    return lowresresult;
		
	}
    else if(selection == "ccurve_in")
	{
	    if(bootstrapiterator == 0)
		{
		    value = highresccurve+".fits";
		    return value;
		}
	    else
		{
		    value = highresccurve+"_BS"+helper4+".fits";
		    return value;
		}
	}

    else if(selection == "msystem_in")
	{
	    if(bootstrapiterator == 0)
		{
		    value = highresmsystemsfile+".fits";
		    return value;
		}
	    else
		{
		    value = highresmsystemsfile+"_BS"+helper4+".fits";
		    return value;
		}
	}

   else if(selection == "msystem_reset")
	{
	    if(bootstrapiterator == 0)
		{
		    value = msystemsoutput+".fits";
		    return value;
		}
	    else
		{
		    value = msystemsoutput+"_BS"+helper4+".fits";
		    return value;
		}
	}

   else if(selection == "core_out")
	{
	    if(bootstrapiterator == 0)
		{
		    value = interpolationoutput+".fits";
		    return value;
		}
	    else
		{
		    value = interpolationoutput+"_BS"+helper4+".fits";
		    return value;
		}
	}

   else if(selection == "merge_out")
	{
	    if(bootstrapiterator == 0)
		{
		    value = mergedoutput+".fits";
		    return value;
		}
	    else
		{
		    value = mergedoutput+"_BS"+helper4+".fits";
		    return value;
		}
	}
   else
       {
	   throw invalid_argument("Invalid selection for show_filename");
       }

}

void HighresOptions::iterate()
{
    bootstrapiterator++;
}




 

BootstrapFieldOptions::BootstrapFieldOptions(const string &filename)
{

  ifstream parafile (filename.c_str());
  string paraline;
  string parameter;

  getline(parafile,paraline);
  istringstream cut1(parametercut(paraline));
  cut1 >>bsstart;
  cut1.str("");

  getline(parafile,paraline);
  istringstream cut2(parametercut(paraline));
  cut2 >>bsstop;
  cut2.str("");

  getline(parafile,paraline);
  istringstream cut3(parametercut(paraline));
  cut3 >>startres;
  cut3.str("");

  getline(parafile,paraline);
  istringstream cut4(parametercut(paraline));
  cut4 >>stopres;
  cut4.str("");

  getline(parafile,paraline);
  istringstream cut5(parametercut(paraline));
  cut5 >>stepsize;
  cut5.str("");

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  square = parameter == "Yes";

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  weight = parameter == "Yes";

  getline(parafile,paraline);
  ellipinput = parametercut(paraline);

  getline(parafile,paraline);
  maskinput = parametercut(paraline);

  getline(parafile,paraline);
  outputdir = parametercut(paraline);

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  cuts1 = gsl_vector_calloc(4);
  cutdoubles(parameter,",",cuts1,4);

  getline(parafile,paraline);
  parameter = parametercut(paraline);
  cuts2 = gsl_vector_calloc(4);
  cutdoubles(parameter,",",cuts2,4);

  getline(parafile,paraline);
  istringstream cut6(parametercut(paraline));
  cut6 >>numgal;
  cut6.str("");

  getline(parafile,paraline);
  istringstream cut7(parametercut(paraline));
  cut7 >>radiusinc;
  cut7.str("");

  getline(parafile,paraline);
  istringstream cut8(parametercut(paraline));
  cut8 >>pixelradius;
  cut8.str("");

  numofsteps = (int) ceil((double)(stopres-startres)/stepsize)+1;

  steps = gsl_vector_int_calloc(numofsteps);
  for(int i = 0; i < (int) ceil((double)(stopres-startres)/stepsize); i++)
    {
      gsl_vector_int_set(steps,i,startres+i*stepsize);
    }
  gsl_vector_int_set(steps,(int) ceil((double)(stopres-startres)/stepsize),stopres);

}

BootstrapFieldOptions::~BootstrapFieldOptions()
{

  gsl_vector_free(cuts1);
  gsl_vector_free(cuts2);
  ellipinput = "";
  maskinput = "";
  outputdir = "";
  gsl_vector_int_free(steps);

}

int BootstrapFieldOptions::showbootstrapindex(int selection)
{

  if(selection == 0)
    {
      return bsstart;
    }
  else
    {
      return bsstop;
    }
}

int BootstrapFieldOptions::showresolution(int selection)
{

  if(selection == 0)
    {
      return startres;
    }
  else
    {
      return stopres;
    }
}

int BootstrapFieldOptions::showstepsize()
{
  return stepsize;
}

int BootstrapFieldOptions::showstepbyindex(int index)
{
  if(index <= numofsteps)
    {
      return gsl_vector_int_get(steps,index);
    }

  else
    {
      return 0;
    }
}

int BootstrapFieldOptions::shownumofsteps()
{
  return numofsteps;
}

bool BootstrapFieldOptions::showflag(int selection)
{
  if(selection == 0)
    {
      return square;
    }
  else
    {
      return weight;
    }
}

string BootstrapFieldOptions::showfilename(int selection)
{
  if(selection == 0)
    {
      return ellipinput;
    }
  else if(selection == 1)
    {
      return outputdir;
    }
  else if(selection == 2)
    {
      return maskinput;
    }
}

double BootstrapFieldOptions::showcuts(int selection)
{
  if(selection == 0)
    {
      return gsl_vector_get(cuts1,0);
    }
  else if(selection == 1)
    {
      return gsl_vector_get(cuts1,1);
    }
  else if(selection == 2)
    {
      return gsl_vector_get(cuts1,2);
    }
  else if(selection == 3)
    {
      return gsl_vector_get(cuts1,3);
    }

  else if(selection == 4)
    {
      return gsl_vector_get(cuts2,0);
    }
  else if(selection == 5)
    {
      return gsl_vector_get(cuts2,1);
    }
  else if(selection == 6)
    {
      return gsl_vector_get(cuts2,2);
    }
  else if(selection == 7)
    {
      return gsl_vector_get(cuts2,3);
    }
}

int BootstrapFieldOptions::showgal()
{
  return numgal;
}

double BootstrapFieldOptions::showradius()
{
  return radiusinc;
}

int BootstrapFieldOptions::showpixelradius()
{
  return pixelradius;
}


BootstrapReconstructionOptions::BootstrapReconstructionOptions(const string &filename)
{

  ifstream parafile (filename.c_str());
  string paraline;
  string parameter;

  getline(parafile,paraline);
  istringstream cut1(parametercut(paraline));
  cut1 >>bsstart;
  cut1.str("");

  getline(parafile,paraline);
  istringstream cut2(parametercut(paraline));
  cut2 >>bsstop;
  cut2.str("");

  getline(parafile,paraline);
  istringstream cut3(parametercut(paraline));
  cut3 >>startdim;
  cut3.str("");

  getline(parafile,paraline);
  istringstream cut4(parametercut(paraline));
  cut4 >>stopdim;
  cut4.str("");

  getline(parafile,paraline);
  istringstream cut5(parametercut(paraline));
  cut5 >>stepsize;
  cut5.str("");

  getline(parafile,paraline);
  infile = parametercut(paraline);

  getline(parafile,paraline);
  outfile = parametercut(paraline);

  getline(parafile,paraline);
  istringstream cut6(parametercut(paraline));
  cut6 >>regparameter;
  cut6.str("");

  getline(parafile,paraline);
  istringstream cut7(parametercut(paraline));
  cut7 >>galredshift;
  cut7.str("");

  getline(parafile,paraline);
  istringstream cut8(parametercut(paraline));
  cut8 >>clusterredshift;
  cut8.str("");

  numofsteps = (int) ceil((double)(stopdim-startdim)/stepsize)+1;

  steps = gsl_vector_int_calloc(numofsteps);
  for(int i = 0; i < (int) ceil((double)(stopdim-startdim)/stepsize); i++)
    {
      gsl_vector_int_set(steps,i,startdim+i*stepsize);
    }
  gsl_vector_int_set(steps,(int) ceil((double)(stopdim-startdim)/stepsize),stopdim);

}


BootstrapReconstructionOptions::~BootstrapReconstructionOptions()
{
  infile = "";
  outfile= "";
  gsl_vector_int_free(steps);
}

int BootstrapReconstructionOptions::showbootstrapindex(int selection)
{
  if(selection == 0)
    {
      return bsstart;
    }
  else
    {
      return bsstop;
    }
}

int BootstrapReconstructionOptions::showresolution(int selection)
{
  if(selection == 0)
    {
      return startdim;
    }
  else
    {
      return stopdim;
    }
}


int BootstrapReconstructionOptions::showstepsize()
{
  return stepsize;
}

int BootstrapReconstructionOptions::showstepbyindex(int index)
{
  return gsl_vector_int_get(steps,index);
}

int BootstrapReconstructionOptions::shownumofsteps()
{
  return numofsteps;
}

string  BootstrapReconstructionOptions::showfilename(int selection)
{
  if(selection == 0)
    {
      return infile;
    }
  else
    {
      return outfile;
    }
}


double BootstrapReconstructionOptions::showreg()
{
  return regparameter;
}


double BootstrapReconstructionOptions::showredshift(int selection)
{
  if(selection == 0)
    {
      return galredshift;
    }

  else
    {
      return clusterredshift;
    }
}

BootstrapAnalysisOptions::BootstrapAnalysisOptions(const string &filename)
{

  ifstream parafile (filename.c_str());
  string paraline;
  string parameter;

  getline(parafile,paraline);
  istringstream cut1(parametercut(paraline));
  cut1 >>realisations;
  cut1.str("");

  getline(parafile,paraline);
  indir = parametercut(paraline);

  getline(parafile,paraline);
  outdir = parametercut(paraline);

  getline(parafile,paraline);
  istringstream cut2(parametercut(paraline));
  cut2 >>initres;
  cut2.str("");

  getline(parafile,paraline);
  istringstream cut3(parametercut(paraline));
  cut3 >>interres;
  cut3.str("");

}


BootstrapAnalysisOptions::~BootstrapAnalysisOptions()
{

  indir = "";
  outdir = "";

}

int BootstrapAnalysisOptions::showrealisations()
{
  return realisations;
}

string BootstrapAnalysisOptions::showfilename(int selection)
{

  if(selection == 0)
    {
      return indir;
    }

  else
    {
      return outdir;
    }
}

int BootstrapAnalysisOptions::showresolution(int selection)
{

  if(selection == 0)
    {
      return initres;
    }

  else
    {
      return interres;
    }

}
 





















 





  
  


























