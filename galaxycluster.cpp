#include "galaxycluster.h"

using namespace std;


GalaxyCluster::GalaxyCluster(gsl_matrix_int *fieldmaskgiven, double x_frctgiven)
{
  fieldmask = gsl_matrix_int_calloc(fieldmaskgiven->size1,fieldmaskgiven->size2);
  gsl_matrix_int_memcpy(fieldmask,fieldmaskgiven);
  x_frct = x_frctgiven;
  FinDifGrid clustergrid(fieldmask,x_frct,true);
  x_dim = clustergrid.showint("x_dim");
  y_dim = clustergrid.showint("y_dim");
  fieldpixels = clustergrid.showint("fieldpixels");

  pot = gsl_vector_calloc(fieldpixels);
  a1 = gsl_vector_calloc(fieldpixels);
  a2 = gsl_vector_calloc(fieldpixels);
  shear1 = gsl_vector_calloc(fieldpixels);
  shear2 = gsl_vector_calloc(fieldpixels);
  convergence = gsl_vector_calloc(fieldpixels);
  f1 = gsl_vector_calloc(fieldpixels);
  f2 = gsl_vector_calloc(fieldpixels);
  g1 = gsl_vector_calloc(fieldpixels);
  g2 = gsl_vector_calloc(fieldpixels);
  jacdet = gsl_vector_calloc(fieldpixels);

  
}

GalaxyCluster::~GalaxyCluster()
{

  gsl_matrix_int_free(fieldmask);
  gsl_vector_free(pot);
  gsl_vector_free(a1);
  gsl_vector_free(a2);
  gsl_vector_free(shear1);
  gsl_vector_free(shear2);
  gsl_vector_free(convergence);
  gsl_vector_free(f1);
  gsl_vector_free(f2);
  gsl_vector_free(g1);
  gsl_vector_free(g2);
  gsl_vector_free(jacdet);

}

int GalaxyCluster::showint(string selection)
{
  if(selection == "x_dim")
    {
      return x_dim;
    }
  else if(selection == "y_dim")
    {
      return y_dim;
    }
  else if(selection == "dim")
    {
      return x_dim*y_dim;
    }
  else if(selection == "fieldpixels")
    {
      return fieldpixels;
    }
  else
    {
      throw invalid_argument("Selection not valid for showint");
    }
}

double GalaxyCluster::showdouble(string selection)
{

  if(selection == "x_frct")
    {
      return x_frct;
    }
  else
    {
      throw invalid_argument("Selection not valid for showdouble");
    }
}

gsl_vector* GalaxyCluster::data(string selection)
{
  if(selection == "pot")
    {
      return pot;
    }
  else if(selection == "shear1")
    {
      return shear1;
    }
  else if(selection == "shear2")
    {
      return shear2;
    }
  else if(selection == "convergence")
    {
      return convergence;
    }
  else if(selection == "f1")
    {
      return f1;
    }
  else if(selection == "f2")
    {
      return f2;
    }
  else if(selection == "g1")
    {
      return g1;
    }
  else if(selection == "g2")
    {
      return g2;
    }
  else if(selection == "jacdet")
    {
      return jacdet;
    }
  else
    {
      throw invalid_argument("Selection not valid for data");
    }
}

gsl_matrix_int* GalaxyCluster::grid()
{
  return fieldmask;
}

void GalaxyCluster::buildfrompot()
{

  //mostly just applying the advanced routines of the FinDifGrid class

  FinDifGrid clustergrid(fieldmask,x_frct,true);
  clustergrid.a1multvec_fast(pot,a1);
  clustergrid.a2multvec_fast(pot,a2);
  clustergrid.s1multvec_fast(pot,shear1);
  clustergrid.s2multvec_fast(pot,shear2);
  clustergrid.cmultvec_fast(pot,convergence);
  clustergrid.f1multvec_fast(pot,f1);
  clustergrid.f2multvec_fast(pot,f2);
  clustergrid.g1multvec_fast(pot,g1);
  clustergrid.g2multvec_fast(pot,g2);

  for(int i = 0; i < fieldpixels; i++)    
    {
      gsl_vector_set(jacdet,i,((1.0-gsl_vector_get(convergence, i))*(1.0-gsl_vector_get(convergence, i)) - (gsl_vector_get(shear1,i)*gsl_vector_get(shear1,i) + gsl_vector_get(shear2,i)*gsl_vector_get(shear2,i))));
    }
 

}

void GalaxyCluster::buildwithoutpot()
{
  for(int i = 0; i < fieldpixels; i++)    
    {
      gsl_vector_set(jacdet,i,((1.0-gsl_vector_get(convergence, i))*(1.0-gsl_vector_get(convergence, i)) - (gsl_vector_get(shear1,i)*gsl_vector_get(shear1,i) + gsl_vector_get(shear2,i)*gsl_vector_get(shear2,i))));
    }
}


void GalaxyCluster::masssheetnormalise(string selection)
{
 
  //normalises the field by setting the pixel with the lowest
  //convergence value to zero

  if(selection == "zero")
    {
      double minimum = gsl_vector_min(convergence);
      double lambda = 1.0/(1.0-minimum);
      gsl_vector_scale(convergence,lambda);
      gsl_vector_add_constant(convergence,1.0-lambda);
      gsl_vector_scale(shear1,lambda);
      gsl_vector_scale(shear2,lambda);
    }
}

void GalaxyCluster::readgridvector(gsl_vector *input,string selection)
{
  if(selection == "pot")
    {
      gsl_vector_memcpy(pot,input);
    }
  else if(selection == "shear1")
    {
      gsl_vector_memcpy(shear1,input);
    }
  else if(selection == "shear2")
    {
      gsl_vector_memcpy(shear2,input);
    }
  else if(selection == "convergence")
    {
      gsl_vector_memcpy(convergence,input);
    }
  else if(selection == "f1")
    {
      gsl_vector_memcpy(f1,input);
    }
  else if(selection == "f2")
    {
      gsl_vector_memcpy(f2,input);
    }
  else if(selection == "g1")
    {
      gsl_vector_memcpy(g1,input);
    }
  else if(selection == "g2")
    {
      gsl_vector_memcpy(g2,input);
    }
  else if(selection == "jacdet")
    {
      gsl_vector_memcpy(jacdet,input);
    }
  else
    {
      throw invalid_argument("Selection not valid for readgrdivector");
    }
}

void GalaxyCluster::readmap(gsl_matrix *input,string selection)
{
  if(selection == "pot")
    {
      maptogridvector(input,pot);
    }
  else if(selection == "shear1")
    {
      maptogridvector(input,shear1);
    }
  else if(selection == "shear2")
    {
      maptogridvector(input,shear2);
    }
  else if(selection == "convergence")
    {
      maptogridvector(input,convergence);
    }
  else if(selection == "f1")
    {
      maptogridvector(input,f1);
    }
  else if(selection == "f2")
    {
      maptogridvector(input,f2);
    }
  else if(selection == "g1")
    {
      maptogridvector(input,g1);
    }
  else if(selection == "g2")
    {
      maptogridvector(input,g2);
    }
  else if(selection == "jacdet")
    {
      maptogridvector(input,jacdet);
    }
  else
    {
      throw invalid_argument("Selection not valid for readmap");
    }
}

void GalaxyCluster::writegridvector(gsl_vector *output, string selection)
{
  if(selection == "pot")
    {
      gsl_vector_memcpy(output,pot);
    }
  else if(selection == "shear1")
    {
      gsl_vector_memcpy(output,shear1);
    }
  else if(selection == "shear2")
    {
      gsl_vector_memcpy(output,shear2);
    }
  else if(selection == "convergence")
    {
      gsl_vector_memcpy(output,convergence);
    }
  else if(selection == "f1")
    {
      gsl_vector_memcpy(output,f1);
    }
  else if(selection == "f2")
    {
      gsl_vector_memcpy(output,f2);
    }
  else if(selection == "g1")
    {
      gsl_vector_memcpy(output,g1);
    }
  else if(selection == "g2")
    {
      gsl_vector_memcpy(output,g2);
    }
  else if(selection == "jacdet")
    {
      gsl_vector_memcpy(output,jacdet);
    }
  else
    {
      throw invalid_argument("Selection not valid for writegridvector");
    }
}

void GalaxyCluster::writemap(gsl_matrix *output, string selection)
{

  // uses another method of this class, described later

  if(selection == "pot")
    {
      gridvectortomap(pot,output);
    }
  else if(selection == "shear1")
    {
      gridvectortomap(shear1,output);
    }
  else if(selection == "shear2")
    {
      gridvectortomap(shear2,output);
    }
  else if(selection == "convergence")
    {
      gridvectortomap(convergence,output);
    }
  else if(selection == "f1")
    {
      gridvectortomap(f1,output);
    }
  else if(selection == "f2")
    {
      gridvectortomap(f2,output);
    }
  else if(selection == "g1")
    {
      gridvectortomap(g1,output);
    }
  else if(selection == "g2")
    {
      gridvectortomap(g2,output);
    }
  else if(selection == "jacdet")
    {
      gridvectortomap(jacdet,output);
    }
  else
    {
      throw invalid_argument("Selection not valid for writemap");
    }
}

void GalaxyCluster::gridvectortomap(gsl_vector *input,gsl_matrix *output)
{
  // Uses again just the routines of the FinDifGrid "working horse"
  FinDifGrid clustergrid(fieldmask,x_frct,true);
  clustergrid.gridvectortomap(input,output,0.0);
}

void GalaxyCluster::maptogridvector(gsl_matrix *input, gsl_vector *output)
{
  // Uses again just the routines of the FinDifGrid "working horse"
  FinDifGrid clustergrid(fieldmask,x_frct,true);
  clustergrid.maptogridvector(input,output);
}

void GalaxyCluster::gridvectortomap(gsl_vector_int *input,gsl_matrix_int *output)
{
  // Uses again just the routines of the FinDifGrid "working horse"
  FinDifGrid clustergrid(fieldmask,x_frct,true);
  clustergrid.gridvectortomap(input,output,100);
}

void GalaxyCluster::maptogridvector(gsl_matrix_int *input, gsl_vector_int *output)
{  
  // Uses again just the routines of the FinDifGrid "working horse"
  FinDifGrid clustergrid(fieldmask,x_frct,true);
  clustergrid.maptogridvector(input,output);
}

void GalaxyCluster::writetofits(string filename)
{
  struct tm *date;
  time_t tim;
  time(&tim);
  
  date = localtime(&tim);
  string t = asctime(date);
  t = t.substr(0,24);

  // the output needs to be in map format for visualisation
  gsl_matrix *potential = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *conv = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *s1 = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *s2 = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *G1 = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *G2 = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *F1 = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *F2 = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *JacDet = gsl_matrix_calloc(y_dim,x_dim);

  //converting the gridvectors to maps
  gridvectortomap(pot,potential);
  gridvectortomap(convergence,conv);
  gridvectortomap(shear1,s1);
  gridvectortomap(shear2,s2);
  gridvectortomap(f1,F1);
  gridvectortomap(f2,F2);
  gridvectortomap(g1,G1);
  gridvectortomap(g2,G2);
  gridvectortomap(jacdet,JacDet);

  //using FITS routines to write to file
  write_pimg(filename,potential);
  write_imge(filename,"convergence",conv);
  write_imge(filename,"shear1",s1);
  write_imge(filename,"shear2",s2);
  write_imge(filename,"f1",G1);
  write_imge(filename,"f2",G2);
  write_imge(filename,"g1",G1);
  write_imge(filename,"g2",G2);
  write_imge(filename,"jacdet",JacDet);
  write_imgeint(filename,"field_mask", fieldmask);

  gsl_matrix_free(potential);
  gsl_matrix_free(conv);
  gsl_matrix_free(s1);
  gsl_matrix_free(s2);
  gsl_matrix_free(G1);
  gsl_matrix_free(G2);
  gsl_matrix_free(F1);
  gsl_matrix_free(F2);
  gsl_matrix_free(JacDet);


}

void GalaxyCluster::writetofits(ReconstructionOptions &options, int iterationindex)
{
  struct tm *date;
  time_t tim;
  time(&tim);
  
  date = localtime(&tim);
  string t = asctime(date);
  t = t.substr(0,24);

  string filename = options.show_filename("output",iterationindex);
  // the output needs to be in map format for visualisation
  gsl_matrix *potential = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *conv = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *s1 = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *s2 = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *G1 = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *G2 = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *F1 = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *F2 = gsl_matrix_calloc(y_dim,x_dim);
  gsl_matrix *JacDet = gsl_matrix_calloc(y_dim,x_dim);

  //converting the gridvectors to maps
  gridvectortomap(pot,potential);
  gridvectortomap(convergence,conv);
  gridvectortomap(shear1,s1);
  gridvectortomap(shear2,s2);
  gridvectortomap(f1,F1);
  gridvectortomap(f2,F2);
  gridvectortomap(g1,G1);
  gridvectortomap(g2,G2);
  gridvectortomap(jacdet,JacDet);

  //using FITS routines to write to file
  write_pimg(filename,potential);
  write_imge(filename,"convergence",conv);
  write_imge(filename,"shear1",s1);
  write_imge(filename,"shear2",s2);
  write_imge(filename,"f1",G1);
  write_imge(filename,"f2",G2);
  write_imge(filename,"g1",G1);
  write_imge(filename,"g2",G2);
  write_imge(filename,"jacdet",JacDet);
  write_imgeint(filename,"field_mask", fieldmask);

  write_header(filename,"DATE",t,"Creation time");
  string method = " ";
  if(options.showbool("shear"))
    {
      method += "s ";
    }
  if(options.showbool("flexion"))
    {
      method += "f ";
    }
  if(options.showbool("ccurve"))
    {
      method += "c ";
    }
  if(options.showbool("msystems"))
    {
      method += "m ";
    }
  write_header(filename,"METHOD",method,"Constraints used for rec.");
  write_header(filename,"ORIGINS",options.showstring("shearinput",iterationindex),"shear Field");
  write_header(filename,"ORIGINF",options.showstring("flexioninput",iterationindex),"flexion Field");
  write_header(filename,"ORIGINC",options.showstring("ccurveinput",iterationindex),"ccurve Field");
  write_header(filename,"ORIGINC",options.showstring("msysteminput",iterationindex),"msystem Field");
  write_header(filename,"REGSHEAR",options.showdouble("regshear"),"shear regularisation");
  write_header(filename,"REGFLEXION",options.showdouble("regflexion"),"flexion regularisation");
  write_header(filename,"REGMETHOD",options.showstring("regtype"),"Regularisation scheme");
  write_header(filename,"CLUSTERZ",options.showdouble("clusterz"),"Redshift of cluster");
  write_header(filename,"SOURCEZ",options.showdouble("shearz"),"Avg source redshift");
  gsl_matrix_free(potential);
  gsl_matrix_free(conv);
  gsl_matrix_free(s1);
  gsl_matrix_free(s2);
  gsl_matrix_free(G1);
  gsl_matrix_free(G2);
  gsl_matrix_free(F1);
  gsl_matrix_free(F2);
  gsl_matrix_free(JacDet);
}


































 



  


	  
            
