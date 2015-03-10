#include <iostream>
#include <fstream>
#include <string>
#include <tclap/CmdLine.h>

using namespace std;

int main(int argc, char* argv[])
{

  TCLAP::CmdLine cmd("Cosmic weight calculation", ' ', "0.1");
  TCLAP::ValueArg<std::string> outputArg("o","output_prefix","Output prefix for config files",true,"","string",cmd);

  cmd.parse( argc, argv );


  string littlehelper = " = \" \"";

  string prefix = outputArg.getValue();
  string out1 = prefix + "_field.conf";
  string out2 = prefix + "_rec.conf";
  string out3 = prefix + "_highres.conf";
  string out4 = prefix + "_mask.conf";

  cout <<endl;
  cout <<"Creating: " <<out1 <<"..." <<flush;
  ofstream outfile1(out1.c_str());

  outfile1 <<"Square field" <<littlehelper <<endl;
  outfile1 <<"Calculate shear" <<littlehelper <<endl;
  outfile1 <<"Calculate flexion" <<littlehelper <<endl;
  outfile1 <<"Calculate ccurve" <<littlehelper <<endl;
  outfile1 <<"Caluclate multiple image positions" <<littlehelper <<endl;
  outfile1 <<"Number of multiple image systems" <<littlehelper <<endl;
  outfile1 <<"Ellipticity input" <<littlehelper <<endl;
  outfile1 <<"Flexion input" <<littlehelper <<endl;
  outfile1 <<"CCurve input" <<littlehelper <<endl;
  outfile1 <<"Multiple image input" <<littlehelper <<endl;
  outfile1 <<"Mask input" <<littlehelper <<endl;
  outfile1 <<"Output prefix for ellip files" <<littlehelper <<endl;
  outfile1 <<"Output prefix for ccurve files" <<littlehelper <<endl;
  outfile1 <<"Output prefix for flexion files" <<littlehelper <<endl;
  outfile1 <<"Output prefix for multiple image files" <<littlehelper <<endl;
  outfile1 <<"Cluster name" <<littlehelper <<endl;
  outfile1 <<"Area of used galaxies" <<littlehelper <<endl;
  outfile1 <<"Area of output field" <<littlehelper <<endl;
  outfile1 <<"Number of galaxies to average for ellip" <<littlehelper <<endl;
  outfile1 <<"Number of galaxies to average for flexion" <<littlehelper <<endl;
  outfile1 <<"Iteration start value" <<littlehelper <<endl;
  outfile1 <<"Iteration stop value" <<littlehelper <<endl;
  outfile1 <<"Stepsize" <<littlehelper <<endl;
  outfile1 <<"Relative radius increment" <<littlehelper <<endl;
  outfile1 <<"Covariance pixelradius" <<littlehelper <<endl;

  outfile1.close();
  cout <<"Done" <<endl;

  cout <<endl;
  cout <<"Creating: " <<out2 <<"..." <<flush;
  ofstream outfile2(out2.c_str());

  outfile2 <<"Use shear" <<littlehelper <<endl;
  outfile2 <<"Use flexion" <<littlehelper <<endl;
  outfile2 <<"Use ccurve" <<littlehelper <<endl;
  outfile2 <<"Use multiple image systems" <<littlehelper <<endl;
  outfile2 <<"Number of multiple image systems" <<littlehelper <<endl;
  outfile2 <<"Iteration start value" <<littlehelper <<endl;
  outfile2 <<"Iteration stop value" <<littlehelper <<endl;
  outfile2 <<"Increment step size" <<littlehelper <<endl;
  outfile2 <<"Set manual steps" <<littlehelper <<endl;
  outfile2 <<"Steps" <<littlehelper <<endl;
  outfile2 <<"Inner level iteration steps" <<littlehelper <<endl;
  outfile2 <<"Convergence threshold for iteration stop" <<littlehelper <<endl;
  outfile2 <<"Prefix of shear-field files" <<littlehelper <<endl;
  outfile2 <<"Prefix of flexion-field files" <<littlehelper <<endl;
  outfile2 <<"Prefix of ccurve files" <<littlehelper <<endl;
  outfile2 <<"Prefix of multiple image system files" <<littlehelper <<endl;
  outfile2 <<"Prefix of output files" <<littlehelper <<endl;
  outfile2 <<"Regularisation type" <<littlehelper <<endl;
  outfile2 <<"Shear regularisation parameter" <<littlehelper <<endl;
  outfile2 <<"Flexion regularisation parameter" <<littlehelper <<endl;
  outfile2 <<"Average redshift of the shear background sources" <<littlehelper <<endl;
  outfile2 <<"Average redshift of the flexion background sources" <<littlehelper <<endl;
  outfile2 <<"Redshift of the cluster" <<littlehelper <<endl;

  outfile2.close();
  cout <<"Done" <<endl;

  cout <<endl;
  cout <<"Creating: " <<out3 <<"..." <<flush;
  ofstream outfile3(out3.c_str());

  outfile3 <<"Use critical curve information" <<littlehelper <<endl;
  outfile3 <<"Use multiple image systems" <<littlehelper <<endl;
  outfile3 <<"Filename of lowres result" <<littlehelper <<endl;
  outfile3 <<"Filename for high res ccurve info" <<littlehelper <<endl;
  outfile3 <<"Filename for high res multiple image system info" <<littlehelper <<endl;
  outfile3 <<"Filename for interpolated clustercore" <<littlehelper <<endl;
  outfile3 <<"Filenname for reset multiple image system info" <<littlehelper <<endl;
  outfile3 <<"Filename for merged maps" <<littlehelper <<endl;
  outfile3 <<"High resolution x-dim" <<littlehelper <<endl;
  outfile3 <<"High resolution y-dim" <<littlehelper <<endl;
  outfile3 <<"Cutpoints" <<littlehelper <<endl;
  outfile3 <<"Number of multiple image systems" <<littlehelper <<endl;
  outfile3 <<"Regularisation parameter" <<littlehelper <<endl;
  outfile3 <<"Regularisation scheme" <<littlehelper <<endl;
  outfile3 <<"Redshift of cluster" <<littlehelper <<endl;
  outfile3 <<"Strong lensing sigma" <<littlehelper <<endl;

  outfile3.close();
  cout <<"Done" <<endl;

  cout <<endl;
  cout <<"Creating: " <<out4 <<"..." <<flush;
  ofstream outfile4(out4.c_str());

  outfile4 <<"Number of rectangular masks" <<littlehelper <<endl;
  for(int i = 1; i <= 20; i++)
    {
      outfile4 <<"Rec-Mask" <<i <<"-Parameters" <<littlehelper <<endl;
    } 

  outfile4.close();
  cout <<"Done" <<endl;


  return 0;

}
