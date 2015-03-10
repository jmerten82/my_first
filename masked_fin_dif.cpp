#include "masked_fin_dif.h"

FinDifGrid::FinDifGrid(gsl_matrix_int *fieldmask, double x_frctgiven, bool switchy)

{

  x_frct = x_frctgiven;
  x_dim = fieldmask->size2;
  y_dim = fieldmask->size1;
  gridmask = gsl_matrix_int_calloc(y_dim,x_dim);
  typemap = gsl_matrix_int_calloc(y_dim,x_dim);
  thirdorder = switchy;
  gsl_matrix_int_memcpy(gridmask,fieldmask);
  badpixels = 0;
  fieldpixels = 0;
  bool badgrid = true;


  //Assigning type to each pixel

  bool notfoundyet;
  if(thirdorder)
    {
      while(badgrid)
	{
	  badpixels = 0;
	  fieldpixels = 0;
	  for(int i = 0; i < y_dim; i++)
	    {
	      for(int j = 0; j < x_dim; j++)
		{
		  notfoundyet = true;
		  if(gsl_matrix_int_get(gridmask,i,j) != 1)
		    {
		      if(notfoundyet && i+2 < y_dim && i-2 >= 0 && j+2 < x_dim && j-2 >= 0)
			{
			  if(gsl_matrix_int_get(gridmask,i+2,j) != 1 && gsl_matrix_int_get(gridmask,i+1,j) !=1 && gsl_matrix_int_get(gridmask,i-1,j) !=1 &&gsl_matrix_int_get(gridmask,i-2,j) != 1 && gsl_matrix_int_get(gridmask,i,j-2) != 1 && gsl_matrix_int_get(gridmask,i,j-1) != 1 && gsl_matrix_int_get(gridmask,i,j+1) != 1 && gsl_matrix_int_get(gridmask,i,j+2) != 1 && gsl_matrix_int_get(gridmask,i+1,j-1) != 1 && gsl_matrix_int_get(gridmask,i+1,j+1) != 1 && gsl_matrix_int_get(gridmask,i-1,j-1) != 1 && gsl_matrix_int_get(gridmask,i-1,j+1) != 1)
			    {
			      gsl_matrix_int_set(typemap,i,j,9);
			      notfoundyet = false;
			      fieldpixels++;
			    } 
			}
		      if(notfoundyet && i+2 < y_dim && i-2 >= 0 && j-3 >= 0)
			{
			  if(gsl_matrix_int_get(gridmask,i+2,j) != 1 && gsl_matrix_int_get(gridmask,i+1,j) !=1 && gsl_matrix_int_get(gridmask,i-1,j) !=1 &&gsl_matrix_int_get(gridmask,i-2,j) != 1 && gsl_matrix_int_get(gridmask,i,j-3) != 1 && gsl_matrix_int_get(gridmask,i,j-2) != 1 && gsl_matrix_int_get(gridmask,i,j-1) != 1 && gsl_matrix_int_get(gridmask,i+1,j-2) != 1 && gsl_matrix_int_get(gridmask,i+1,j-1) != 1 && gsl_matrix_int_get(gridmask,i-1,j-2) != 1 && gsl_matrix_int_get(gridmask,i-1,j-1) != 1)
			    {
			      gsl_matrix_int_set(typemap,i,j,8);
			      notfoundyet = false;
			      fieldpixels++;
			    }
			}
		      if(notfoundyet && i+2 < y_dim && i-2 >= 0 && j+3 < x_dim)
			{
			  if(gsl_matrix_int_get(gridmask,i+2,j) != 1 && gsl_matrix_int_get(gridmask,i+1,j) !=1 && gsl_matrix_int_get(gridmask,i-1,j) !=1 &&gsl_matrix_int_get(gridmask,i-2,j) != 1 && gsl_matrix_int_get(gridmask,i,j+3) != 1 && gsl_matrix_int_get(gridmask,i,j+2) != 1 && gsl_matrix_int_get(gridmask,i,j+1) != 1 && gsl_matrix_int_get(gridmask,i+1,j+2) != 1 && gsl_matrix_int_get(gridmask,i+1,j+1) != 1 && gsl_matrix_int_get(gridmask,i-1,j+2) != 1 && gsl_matrix_int_get(gridmask,i-1,j+1) != 1)
			    {
			      gsl_matrix_int_set(typemap,i,j,7);
			      notfoundyet = false;
			      fieldpixels++;
			    }
			}
		      if(notfoundyet && j+2 < x_dim && j-2 >= 0 && i-3 >= 0)
			{
			  if(gsl_matrix_int_get(gridmask,i,j+2) != 1 && gsl_matrix_int_get(gridmask,i,j+1) !=1 && gsl_matrix_int_get(gridmask,i,j-1) !=1 &&gsl_matrix_int_get(gridmask,i,j-2) != 1 && gsl_matrix_int_get(gridmask,i-3,j) != 1 && gsl_matrix_int_get(gridmask,i-2,j) != 1 && gsl_matrix_int_get(gridmask,i-2,j-1) != 1 && gsl_matrix_int_get(gridmask,i-2,j+1) != 1 && gsl_matrix_int_get(gridmask,i-1,j+1) != 1 && gsl_matrix_int_get(gridmask,i-1,j-1) != 1 && gsl_matrix_int_get(gridmask,i-1,j) != 1)
			    {
			      gsl_matrix_int_set(typemap,i,j,6);
			      notfoundyet = false;
			      fieldpixels++;
			    }
			}
		      if(notfoundyet && j+2 < x_dim && j-2 >= 0 && i+3 < y_dim)
			{		    
			  if(gsl_matrix_int_get(gridmask,i,j+2) != 1 && gsl_matrix_int_get(gridmask,i,j+1) !=1 && gsl_matrix_int_get(gridmask,i,j-1) !=1 &&gsl_matrix_int_get(gridmask,i,j-2) != 1 && gsl_matrix_int_get(gridmask,i+3,j) != 1 && gsl_matrix_int_get(gridmask,i+2,j) != 1 && gsl_matrix_int_get(gridmask,i+2,j-1) != 1 && gsl_matrix_int_get(gridmask,i+2,j+1) != 1 && gsl_matrix_int_get(gridmask,i+1,j+1) != 1 && gsl_matrix_int_get(gridmask,i+1,j-1) != 1 && gsl_matrix_int_get(gridmask,i+1,j)!=1)
			    {
			      gsl_matrix_int_set(typemap,i,j,5);
			      notfoundyet = false;
			      fieldpixels++;
			    }		    
			}
		      if(notfoundyet && j-3 >= 0 && i-3 >= 0)
			{
			  if(gsl_matrix_int_get(gridmask,i-3,j) != 1 && gsl_matrix_int_get(gridmask,i-2,j) !=1 && gsl_matrix_int_get(gridmask,i-1,j) !=1 &&gsl_matrix_int_get(gridmask,i-2,j-1) != 1 && gsl_matrix_int_get(gridmask,i-1,j-1) != 1 && gsl_matrix_int_get(gridmask,i,j-1) != 1 && gsl_matrix_int_get(gridmask,i-2,j-2) != 1 && gsl_matrix_int_get(gridmask,i-1,j-2) != 1 && gsl_matrix_int_get(gridmask,i,j-2) != 1 && gsl_matrix_int_get(gridmask,i,j-3) != 1)
			    {
			      gsl_matrix_int_set(typemap,i,j,4);
			      notfoundyet = false;
			      fieldpixels++;
			    }
			}
		      if(notfoundyet && j+3 < x_dim && i-3 >= 0)
			{
			  if(gsl_matrix_int_get(gridmask,i-3,j) != 1 && gsl_matrix_int_get(gridmask,i-2,j) !=1 && gsl_matrix_int_get(gridmask,i-1,j) !=1 &&gsl_matrix_int_get(gridmask,i-2,j+1) != 1 && gsl_matrix_int_get(gridmask,i-1,j+1) != 1 && gsl_matrix_int_get(gridmask,i,j+1) != 1 && gsl_matrix_int_get(gridmask,i-2,j+2) != 1 && gsl_matrix_int_get(gridmask,i-1,j+2) != 1 && gsl_matrix_int_get(gridmask,i,j+2) != 1 && gsl_matrix_int_get(gridmask,i,j+3) != 1)
			    {
			      gsl_matrix_int_set(typemap,i,j,3);
			      notfoundyet = false;
			      fieldpixels++;
			    }
			}
		      
		      if(notfoundyet && j-3 >= 0 && i+3 < y_dim)
			{
			  if(gsl_matrix_int_get(gridmask,i+3,j) != 1 && gsl_matrix_int_get(gridmask,i+2,j) !=1 && gsl_matrix_int_get(gridmask,i+1,j) !=1 &&gsl_matrix_int_get(gridmask,i+2,j-1) != 1 && gsl_matrix_int_get(gridmask,i+1,j-1) != 1 && gsl_matrix_int_get(gridmask,i,j-1) != 1 && gsl_matrix_int_get(gridmask,i+2,j-2) != 1 && gsl_matrix_int_get(gridmask,i+1,j-2) != 1 && gsl_matrix_int_get(gridmask,i,j-2) != 1 && gsl_matrix_int_get(gridmask,i,j-3) != 1)
			    {
			      gsl_matrix_int_set(typemap,i,j,2);
			      notfoundyet = false;
			      fieldpixels++;
			    }
			}
		      if(notfoundyet && j+3 < x_dim && i+3 < y_dim)
			{
			  if(gsl_matrix_int_get(gridmask,i+3,j) != 1 && gsl_matrix_int_get(gridmask,i+2,j) !=1 && gsl_matrix_int_get(gridmask,i+1,j) !=1 &&gsl_matrix_int_get(gridmask,i+2,j+1) != 1 && gsl_matrix_int_get(gridmask,i+1,j+1) != 1 && gsl_matrix_int_get(gridmask,i,j+1) != 1 && gsl_matrix_int_get(gridmask,i+2,j+2) != 1 && gsl_matrix_int_get(gridmask,i+1,j+2) != 1 && gsl_matrix_int_get(gridmask,i,j+2) != 1 && gsl_matrix_int_get(gridmask,i,j+3) != 1)
			    {
			      gsl_matrix_int_set(typemap,i,j,1);
			      notfoundyet = false;
			      fieldpixels++;
			    }
			}
		      if(notfoundyet)
			{
			  badpixels++;
			  gsl_matrix_int_set(typemap,i,j,-1);
			  gsl_matrix_int_set(gridmask,i,j,1);
			}
		    }
		  
		}
	    }
	  if(badpixels == 0)
	    {
	      badgrid = false;
	    }
	}
    }
    
  else
    {
      while(badgrid)
	{
	  badpixels = 0;
	  fieldpixels = 0;
	  for(int i = 0; i < y_dim; i++)
	    {
	      for(int j = 0; j < x_dim; j++)
		{
		  notfoundyet = true;
		  if(gsl_matrix_int_get(gridmask,i,j) != 1)
		    {
		      if(notfoundyet && i+1 < y_dim && i-1 >= 0 && j+1 < x_dim && j-1 >= 0)
			{
			  if(gsl_matrix_int_get(gridmask,i+1,j) !=1 && gsl_matrix_int_get(gridmask,i-1,j) !=1 && gsl_matrix_int_get(gridmask,i,j-1) != 1 && gsl_matrix_int_get(gridmask,i,j+1) != 1 && gsl_matrix_int_get(gridmask,i+1,j-1) != 1 && gsl_matrix_int_get(gridmask,i+1,j+1) != 1 && gsl_matrix_int_get(gridmask,i-1,j-1) != 1 && gsl_matrix_int_get(gridmask,i-1,j+1) != 1)
			    {
			      gsl_matrix_int_set(typemap,i,j,9);
			      notfoundyet = false;
			      fieldpixels++;
			    } 
			} 
		      if(notfoundyet && i+1 < y_dim && i-1 >= 0 && j-2 >= 0)
			{
			  if(gsl_matrix_int_get(gridmask,i+1,j) !=1 && gsl_matrix_int_get(gridmask,i-1,j) !=1 && gsl_matrix_int_get(gridmask,i,j-2) != 1 && gsl_matrix_int_get(gridmask,i,j-1) != 1 && gsl_matrix_int_get(gridmask,i+1,j-1) != 1 && gsl_matrix_int_get(gridmask,i-1,j-1)!=1)
			    {
			      gsl_matrix_int_set(typemap,i,j,8);
			      notfoundyet = false;
			      fieldpixels++;
			    }
			}
		      if(notfoundyet && i+1 < y_dim && i-1 >= 0 && j+2 < x_dim)
			{
			  if(gsl_matrix_int_get(gridmask,i+1,j) !=1 && gsl_matrix_int_get(gridmask,i-1,j) !=1 && gsl_matrix_int_get(gridmask,i,j+2) != 1 && gsl_matrix_int_get(gridmask,i,j+1) != 1 && gsl_matrix_int_get(gridmask,i+1,j+1) != 1 && gsl_matrix_int_get(gridmask,i-1,j+1)!=1)
			    {
			      gsl_matrix_int_set(typemap,i,j,7);
			      notfoundyet = false;
			      fieldpixels++;
			    }
			}
		      if(notfoundyet && j+1 < x_dim && j-1 >= 0 && i-2 >= 0)
			{
			  if(gsl_matrix_int_get(gridmask,i,j+1) !=1 && gsl_matrix_int_get(gridmask,i,j-1) !=1 && gsl_matrix_int_get(gridmask,i-2,j) != 1 && gsl_matrix_int_get(gridmask,i-1,j) != 1 && gsl_matrix_int_get(gridmask,i-1,j+1) != 1 && gsl_matrix_int_get(gridmask,i-1,j-1)!=1)
			    {
			      gsl_matrix_int_set(typemap,i,j,6);
			      notfoundyet = false;
			      fieldpixels++;
			    }
			}
		      if(notfoundyet && j+1 < x_dim && j-1 >= 0 && i+2 < y_dim)
			{
			  if(gsl_matrix_int_get(gridmask,i,j+1) !=1 && gsl_matrix_int_get(gridmask,i,j-1) !=1 && gsl_matrix_int_get(gridmask,i+2,j) != 1 && gsl_matrix_int_get(gridmask,i+1,j) != 1 && gsl_matrix_int_get(gridmask,i+1,j+1) != 1 && gsl_matrix_int_get(gridmask,i+1,j-1)!=1)
			    {
			      gsl_matrix_int_set(typemap,i,j,5);
			      notfoundyet = false;
			      fieldpixels++;
			    }
			}
		      if(notfoundyet && j-2 >= 0 && i-2 >= 0)
			{
			  if(gsl_matrix_int_get(gridmask,i,j-2) !=1 && gsl_matrix_int_get(gridmask,i,j-1) !=1 && gsl_matrix_int_get(gridmask,i-1,j-1) != 1 && gsl_matrix_int_get(gridmask,i-1,j) != 1 && gsl_matrix_int_get(gridmask,i-2,j) != 1)
			    {
			      gsl_matrix_int_set(typemap,i,j,4);
			      notfoundyet = false;
			      fieldpixels++;
			    }
			}
		      if(notfoundyet && j+2 < y_dim && i-2 >= 0)
			{
			  if(gsl_matrix_int_get(gridmask,i,j+2) !=1 && gsl_matrix_int_get(gridmask,i,j+1) !=1 && gsl_matrix_int_get(gridmask,i-1,j+1) != 1 && gsl_matrix_int_get(gridmask,i-1,j) != 1 && gsl_matrix_int_get(gridmask,i-2,j) != 1)
			    {
			      gsl_matrix_int_set(typemap,i,j,3);
			      notfoundyet = false;
			      fieldpixels++;
			    }
			}
		      
		      if(notfoundyet && j-2 >= 0 && i+2 < y_dim)
			{
			  if(gsl_matrix_int_get(gridmask,i+2,j) !=1 && gsl_matrix_int_get(gridmask,i+1,j) !=1 && gsl_matrix_int_get(gridmask,i+1,j-1) != 1 && gsl_matrix_int_get(gridmask,i,j-2) != 1 && gsl_matrix_int_get(gridmask,i,j-1) != 1)
			    {
			      gsl_matrix_int_set(typemap,i,j,2);
			      notfoundyet = false;
			      fieldpixels++;
			    }
			}
		      if(notfoundyet && j+2 < x_dim && i+2 < y_dim)
			{
			  if(gsl_matrix_int_get(gridmask,i+2,j) !=1 && gsl_matrix_int_get(gridmask,i+1,j) !=1 && gsl_matrix_int_get(gridmask,i+1,j+1) != 1 && gsl_matrix_int_get(gridmask,i,j+2) != 1 && gsl_matrix_int_get(gridmask,i,j+1) != 1)
			    {
			      gsl_matrix_int_set(typemap,i,j,1);
			      notfoundyet = false;
			      fieldpixels++;
			    }
			}
		      if(notfoundyet)
			{
			  badpixels++;
			  gsl_matrix_int_set(typemap,i,j,-1);
			  gsl_matrix_int_set(gridmask,i,j,1);
			}
		    }
		}
	    }
	  if(badpixels == 0)
	    {
	      badgrid = false;
	    }
	}
    }

  typevector = gsl_vector_int_calloc(fieldpixels);
  distancetable = gsl_matrix_int_calloc(fieldpixels,10);

  //Assigning indices to the pixels which won't be masked later

  int counter = 0;
  badpixels = 0;

  for(int i = 0; i < y_dim; i++)
    {
      for(int j = 0; j < x_dim; j++)
	{
	  if(gsl_matrix_int_get(typemap,i,j) > 0)
	    {
	      gsl_matrix_int_set(gridmask,i,j,counter);
	      gsl_vector_int_set(typevector,counter,gsl_matrix_int_get(typemap,i,j));
	      counter++;
	    }
	  else
	    {
	      badpixels++;
	      gsl_matrix_int_set(gridmask,i,j,-1);
	    }
	}
    }

  //Checking the distances to the needed neighbours

  counter = 0;

  for(int i = 0; i < y_dim; i++)
    {
      for(int j = 0; j < x_dim; j++)
	{
	  if(gsl_matrix_int_get(typemap,i,j) > 0)
	    {
	      gsl_matrix_int_set(distancetable,counter,0,counter);
	      gsl_matrix_int_set(distancetable,counter,1,j);
	      gsl_matrix_int_set(distancetable,counter,2,i);
	      gsl_matrix_int_set(distancetable,counter,3,gsl_matrix_int_get(typemap,i,j));
	      if(i+3 < y_dim)
		{
		  if(gsl_matrix_int_get(typemap,i+3,j) > 0)
		    {
		      gsl_matrix_int_set(distancetable,counter,4,gsl_matrix_int_get(gridmask,i+3,j)-gsl_matrix_int_get(gridmask,i,j));		    
		    }
		}
	      if(i+2 < y_dim)
		{
		  if(gsl_matrix_int_get(typemap,i+2,j) > 0)
		    {
		      gsl_matrix_int_set(distancetable,counter,5,gsl_matrix_int_get(gridmask,i+2,j)-gsl_matrix_int_get(gridmask,i,j));		    
		    }
		}
	      if(i+1 < y_dim)
		{
		  if(gsl_matrix_int_get(typemap,i+1,j) > 0)
		    {
		      gsl_matrix_int_set(distancetable,counter,6,gsl_matrix_int_get(gridmask,i+1,j)-gsl_matrix_int_get(gridmask,i,j));		    
		    }
		}
	      if(i-1 >= 0)
		{
		  if(gsl_matrix_int_get(typemap,i-1,j) > 0)
		    {
		      gsl_matrix_int_set(distancetable,counter,7,gsl_matrix_int_get(gridmask,i-1,j)-gsl_matrix_int_get(gridmask,i,j));		    
		    }
		}
	      if(i-2 >= 0)
		{
		  if(gsl_matrix_int_get(typemap,i-2,j) > 0)
		    {
		      gsl_matrix_int_set(distancetable,counter,8,gsl_matrix_int_get(gridmask,i-2,j)-gsl_matrix_int_get(gridmask,i,j));		    
		    }
		}
	      if(i-3 >= 0)
		{
		  if(gsl_matrix_int_get(typemap,i-3,j) > 0)
		    {
		      gsl_matrix_int_set(distancetable,counter,9,gsl_matrix_int_get(gridmask,i-3,j)-gsl_matrix_int_get(gridmask,i,j));		    
		    }
		}
	      counter++;
	    }
	}
    }
}

FinDifGrid::~FinDifGrid()
{
  gsl_matrix_int_free(gridmask);
  gsl_matrix_int_free(typemap);
  gsl_matrix_int_free(distancetable);
  gsl_vector_int_free(typevector);
}

int FinDifGrid::showint(string selection)
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
      return x_dim * y_dim;
    }
  else if(selection == "fieldpixels")
    {
      return fieldpixels;
    }
  else if(selection == "badpixels")
    {
      return badpixels;
    }
  else
    {
      throw invalid_argument("Selection not valid for showint");
    }

}

double FinDifGrid::showdouble(string selection)
{

  if(selection == "x_Frct")
    {

      return x_frct;
    }
  else
    {
      throw invalid_argument("Selection not valid for showdouble");
    }

}

void FinDifGrid::givegridmask(gsl_matrix_int *output)
{

  gsl_matrix_int_memcpy(output,gridmask);

}

void FinDifGrid::givetypemap(gsl_matrix_int *output)
{


  gsl_matrix_int_memcpy(output,typemap);
}

void FinDifGrid::givetypemap(gsl_vector_int *output)
{

  gsl_vector_int_memcpy(output,typevector);

}

void FinDifGrid::givedistancetable(gsl_matrix_int *output)
{

  gsl_matrix_int_memcpy(output,distancetable);

}


void FinDifGrid::givedistancetable(const string &output)
{

  ofstream out(output.c_str());

  for(int i = 0; i < fieldpixels; i++)
    {
      out <<gsl_matrix_int_get(distancetable,i,0) <<"\t" <<gsl_matrix_int_get(distancetable,i,1) <<"\t" <<gsl_matrix_int_get(distancetable,i,2) <<"\t" <<gsl_matrix_int_get(distancetable,i,3) <<"\t" <<gsl_matrix_int_get(distancetable,i,4) <<"\t" <<gsl_matrix_int_get(distancetable,i,5) <<"\t" <<gsl_matrix_int_get(distancetable,i,6) <<"\t" <<gsl_matrix_int_get(distancetable,i,7) <<"\t" <<gsl_matrix_int_get(distancetable,i,8) <<"\t" <<gsl_matrix_int_get(distancetable,i,9) <<endl;
    }

  out.close();
}

void FinDifGrid::maptogridvector(gsl_matrix *input, gsl_vector *output)
{
  int index = 0;

  for(int i = 0; i < y_dim; i++)
    {
      for(int j = 0; j < x_dim; j++)
	{
	  if(gsl_matrix_int_get(typemap,i,j) > 0)
	    {
	      gsl_vector_set(output,index,gsl_matrix_get(input,i,j));
	      index++;
	    }
	}
    }
}

void FinDifGrid::maptogridvector(gsl_matrix_int *input, gsl_vector_int *output)
{
  int index = 0;

  for(int i = 0; i < y_dim; i++)
    {
      for(int j = 0; j < x_dim; j++)
	{
	  if(gsl_matrix_int_get(typemap,i,j) > 0)
	    {
	      gsl_vector_int_set(output,index,gsl_matrix_int_get(input,i,j));
	      index++;
	    }
	}
    }
}

void FinDifGrid::gridvectortomap(gsl_vector *input, gsl_matrix *output, double value)
{
  int index = 0;

  for(int i = 0; i < y_dim; i++)
    {
      for(int j = 0; j < x_dim; j++)
	{
	  if(gsl_matrix_int_get(typemap,i,j) > 0)
	    {
	      gsl_matrix_set(output,i,j,gsl_vector_get(input,index));
	      index++;
	    }
	  else
	    {
	      gsl_matrix_set(output,i,j,value);
	    }
	}
    }
}

void FinDifGrid::gridvectortomap(gsl_vector_int *input, gsl_matrix_int *output, int value)
{
  int index = 0;

  for(int i = 0; i < y_dim; i++)
    {
      for(int j = 0; j < x_dim; j++)
	{
	  if(gsl_matrix_int_get(typemap,i,j) > 0)
	    {
	      gsl_matrix_int_set(output,i,j,gsl_vector_int_get(input,index));
	      index++;
	    }
	  else
	    {
	      gsl_matrix_int_set(output,i,j,value);
	    }
	}
    }
}

double FinDifGrid::a1value(int row, int col)
{
  double x = x_dim/x_frct;
  
  if(gsl_vector_int_get(typevector,row) == 1 || gsl_vector_int_get(typevector,row) == 3 || gsl_vector_int_get(typevector,row) == 7)
    {
      if(col == row)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+1)
	{
	  return 2.0*x;
	}
      else if(col == row+2)
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 2 || gsl_vector_int_get(typevector,row) == 4 || gsl_vector_int_get(typevector,row) == 8)
    {
      if(col == row)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row-1)
	{
	  return -2.0*x;
	}
      else if(col == row-2)
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else
    {
      if(col == row-1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+1)
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
}

double FinDifGrid::a2value(int row, int col)
{

  double x = x_dim/x_frct;

  if(gsl_vector_int_get(typevector,row) == 1 || gsl_vector_int_get(typevector,row) == 2 || gsl_vector_int_get(typevector,row) == 5)
    {
      if(col == row)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row + gsl_matrix_int_get(distancetable,row,6))
	{
	  return 2.0*x;
	}
      else if(col == row + gsl_matrix_int_get(distancetable,row,5))
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 3 || gsl_vector_int_get(typevector,row) == 4 || gsl_vector_int_get(typevector,row) == 6)
    {

	if(col == row)
	  {
	    return 3.0/2.0*x;
	  }
	else if(col == row + gsl_matrix_int_get(distancetable,row,7))
	  {
	    return -2.0*x;
	  }
	else if(col == row + gsl_matrix_int_get(distancetable,row,8))
	  {
	    return 1.0/2.0*x;
	  }
	else
	  {
	    return 0.0;
	  }
    }
  else
    {
      if(col == row + gsl_matrix_int_get(distancetable,row,6))
	{
	  return 1.0/2.0*x;
	}
      else if(col == row + gsl_matrix_int_get(distancetable,row,7))
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
}

double FinDifGrid::s1value(int row, int col)
{
  double x = (x_dim * x_dim)/(x_frct*x_frct);

  if(gsl_vector_int_get(typevector,row) == 1)
    {
      if(col == row+1)
	{
	  return -x;
	}
      else if(col == row+2)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;

	}
    }
  else if(gsl_vector_int_get(typevector,row) == 2)
    {
      if(col == row-1)
	{
	  return -x;
	}
      else if(col == row-2)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return x;
	}
      else if(col == row + gsl_matrix_int_get(distancetable,row,5))
	{
	  return -1.0/2.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 3)
    {
      if(col == row+1)
	{
	  return -x;
	}
      else if(col == row+2)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 4)
    {
      if(col == row-1)
	{
	  return -x;
	}
      else if(col == row-2)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 5)
    {
      if(col == row)
	{
	  return -3.0/2.0*x;
	}
      else if(col ==row-1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 6)
    {
      if(col == row)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row-1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 7)
    {
      if(col ==row)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+1)
	{
	  return -x;
	}
      else if(col == row+2)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -1.0/2.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row)== 8)
    {
      if(col == row)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row-1)
	{
	  return -x;
	}
      else if(col == row-2)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else
    {
      if(col == row-1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
}

double FinDifGrid::s2value(int row, int col)
{
  double x = (x_dim * x_dim)/(x_frct*x_frct);

  if(gsl_vector_int_get(typevector,row) == 1)
    {
      if(col == row)
	{
	  return x;
	}
      else if(col == row+1)
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return x;
	}
      else
	{
	  return 0.0;

	}
    }
  else if(gsl_vector_int_get(typevector,row) == 2)
    {
      if(col == row)
	{
	  return -x;
	}
      else if(col == row-1)
	{
	  return x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return x;
	}
      else if(col == row + gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return -x;
	}
      else 
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 3)
    {
      if(col == row)
	{
	  return -x;
	}
      else if(col == row+1)
	{
	  return x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return -x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 4)
    {
      if(col == row)
	{
	  return x;
	}
      else if(col == row-1)
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 5)
    {
      if(col == row-1)
	{
	  return 1.0/2.0*x;
	}
      else if(col ==row+1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 6)
    {
      if(col == row-1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 7)
    {
      if(col ==row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return -1.0/2.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row)== 8)
    {
      if(col ==row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return 1.0/2.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }
  else
    {
      if(col == row)
	{
	  return x;
	}
      else if(col == row-1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row + gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return 1.0/2.0*x;
	}    
      else
	{
	  return 0.0;
	}
    }
}

double FinDifGrid::cvalue(int row, int col)
{
  double x = (x_dim * x_dim)/(x_frct*x_frct);

  if(gsl_vector_int_get(typevector,row) == 1)
    {
      if(col == row)
	{
	  return x;
	}
      else if(col == row+1)
	{
	  return -x;
	}
      else if(col == row+2)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;

	}
    }
  else if(gsl_vector_int_get(typevector,row) == 2)
    {
      if(col == row)
	{
	  return x;
	}
      else if(col == row-1)
	{
	  return -x;
	}
      else if(col == row-2)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -x;
	}
      else if(col == row + gsl_matrix_int_get(distancetable,row,5))
	{
	  return 1.0/2.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 3)
    {
      if(col == row)
	{
	  return x;
	}

      else if(col == row+1)
	{
	  return -x;
	}
      else if(col == row+2)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 4)
    {
      if(col == row)
	{
	  return x;
	}
      else if(col == row-1)
	{
	  return -x;
	}
      else if(col == row-2)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 5)
    {
      if(col == row)
	{
	  return -1.0/2.0*x;
	}
      else if(col ==row-1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 6)
    {
      if(col == row)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row-1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 7)
    {
      if(col ==row)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+1)
	{
	  return -x;
	}
      else if(col == row+2)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return 1.0/2.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row)== 8)
    {
      if(col == row)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row-1)
	{
	  return -x;
	}
      else if(col == row-2)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else
    {
      if(col == row)
	{
	  return -2.0/3.0*x;
	}
      else if(col == row-1)
	{
	  return -1.0/6.0*x;
	}
      else if(col == row+1)
	{
	  return -1.0/6.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -1.0/6.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -1.0/6.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return 1.0/3.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return 1.0/3.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return 1.0/3.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return 1.0/3.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
}

double FinDifGrid::f1value(int row, int col)
{
  double x = (x_dim * x_dim * x_dim)/(x_frct*x_frct*x_frct);

  if(gsl_vector_int_get(typevector,row) == 1)
    {
      if(col == row)
	{
	  return -x;
	}
      else if(col == row+1)
	{
	  return 2.0*x;
	}
      else if(col == row+2)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+3)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5)+1)
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;

	}
    }
  else if(gsl_vector_int_get(typevector,row) == 2)
    {
      if(col == row)
	{
	  return x;
	}
      else if(col == row-1)
	{
	  return -2.0*x;
	}
      else if(col == row-2)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row-3)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return x;
	}
      else if(col == row + gsl_matrix_int_get(distancetable,row,5))
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5)-1)
	{
	  return -1.0/2.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 3)
    {
      if(col == row)
	{
	  return -x;
	}
      else if(col == row+1)
	{
	  return 2.0*x;
	}
      else if(col == row+2)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+3)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8)+1)
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;

	}
    }
  else if(gsl_vector_int_get(typevector,row) == 4)
    {
      if(col == row)
	{
	  return x;
	}
      else if(col == row-1)
	{
	  return -2.0*x;
	}
      else if(col == row-2)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row-3)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return x;
	}
      else if(col == row + gsl_matrix_int_get(distancetable,row,8))
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8)-1)
	{
	  return -1.0/2.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }

  else if(gsl_vector_int_get(typevector,row) == 5)
    {
      if(col ==row-1)
	{
	  return 1.0/4.0*x;
	}
      else if(col == row-2)
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+1)
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+2)
	{
	  return 1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5)-1)
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5)+1)
	{
	  return 1.0/4.0*x;
	}
      else
	{
	  return 0.0;
	}
    }

  else if(gsl_vector_int_get(typevector,row) == 6)
    {
      if(col ==row-1)
	{
	  return 1.0/4.0*x;
	}
      else if(col == row-2)
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+1)
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+2)
	{
	  return 1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8)-1)
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8)+1)
	{
	  return 1.0/4.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 7)
    {
      if(col ==row)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+2)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+3)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return 1.0/2.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 8)
    {
      if(col ==row)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row-1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row-2)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row-3)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return -1.0/2.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }
  else
    {
      if(col == row-1)
	{
	  return x;
	}
      else if(col == row-2)
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+1)
	{
	  return -x;
	}
      else if(col == row+2)
	{
	  return 1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return 1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return 1.0/4.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
}

double FinDifGrid::f2value(int row, int col)
{
  double x = (x_dim * x_dim * x_dim)/(x_frct*x_frct*x_frct);

  if(gsl_vector_int_get(typevector,row) == 1)
    {
      if(col == row)
	{
	  return -x;
	}
      else if(col == row+1)
	{
	  return x;
	}
      else if(col == row+2)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return 2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+2)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,4))
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;

	}
    }
  else if(gsl_vector_int_get(typevector,row) == 2)
    {
      if(col == row)
	{
	  return -x;
	}
      else if(col == row-1)
	{
	  return x;
	}
      else if(col == row-2)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return 2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-2)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,4))
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;

	}
    }
  else if(gsl_vector_int_get(typevector,row) == 3)
    {
      if(col == row)
	{
	  return x;
	}
      else if(col == row+1)
	{
	  return -x;
	}
      else if(col == row+2)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+2)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,9))
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;

	}
    }
  else if(gsl_vector_int_get(typevector,row) == 4)
    {
      if(col == row)
	{
	  return x;
	}
      else if(col == row-1)
	{
	  return -x;
	}
      else if(col == row-2)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-2)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,9))
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;

	}
    }


  else if(gsl_vector_int_get(typevector,row) == 5)
    {
      if(col == row)
	{
	  return 1.0/2.0*x;
	}
      else if(col ==row-1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row +1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,4))
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 6)
    {
      if(col == row)
	{
	  return -1.0/2.0*x;
	}
      else if(col ==row-1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row +1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,9))
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }

  else if(gsl_vector_int_get(typevector,row) == 7)
    {

      if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+2)
	{
	  return 1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return 1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return 1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+2)
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return -1.0/4.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 8)
    {

      if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-2)
	{
	  return 1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return 1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return 1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-2)
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return -1.0/4.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }

  else
    {
      if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return 1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return 1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return 1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return -1.0/4.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
}

double FinDifGrid::g1value(int row, int col)
{
  double x = (x_dim * x_dim * x_dim)/(x_frct*x_frct*x_frct);

  if(gsl_vector_int_get(typevector,row) == 1)
    {
      if(col == row)
	{
	  return x;
	}
      else if(col == row+2)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+3)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -3.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return 3.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5)+1)
	{
	  return -3.0/2.0*x;
	}
      else
	{
	  return 0.0;

	}
    }
  else if(gsl_vector_int_get(typevector,row) == 2)
    {
      if(col == row)
	{
	  return -x;
	}
      else if(col == row-2)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row-3)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return 3.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return -3.0*x;
	}
      else if(col == row + gsl_matrix_int_get(distancetable,row,5))
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5)-1)
	{
	  return 3.0/2.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 3)
    {
      if(col == row)
	{
	  return x;
	}
      else if(col == row+2)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+3)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -3.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return 3.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8)+1)
	{
	  return -3.0/2.0*x;
	}
      else
	{
	  return 0.0;

	}
    }
  else if(gsl_vector_int_get(typevector,row) == 4)
    {
      if(col == row)
	{
	  return -x;
	}
      else if(col == row-2)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row-3)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return 3.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return -3.0*x;
	}
      else if(col == row + gsl_matrix_int_get(distancetable,row,8))
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8)-1)
	{
	  return 3.0/2.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }

  else if(gsl_vector_int_get(typevector,row) == 5)
    {
      if(col ==row-1)
	{
	  return 5.0/4.0*x;
	}
      else if(col == row-2)
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+1)
	{
	  return -5.0/4.0*x;
	}
      else if(col == row+2)
	{
	  return 1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5)-1)
	{
	  return 3.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5)+1)
	{
	  return -3.0/4.0*x;
	}
      else
	{
	  return 0.0;
	}
    }

  else if(gsl_vector_int_get(typevector,row) == 6)
    {
      if(col ==row-1)
	{
	  return 5.0/4.0*x;
	}
      else if(col == row-2)
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+1)
	{
	  return -5.0/4.0*x;
	}
      else if(col == row+2)
	{
	  return 1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8)-1)
	{
	  return 3.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8)+1)
	{
	  return -3.0/4.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 7)
    {
      if(col ==row)
	{
	  return -7.0/2.0*x;
	}
      else if(col == row+1)
	{
	  return 9.0/2.0*x;
	}
      else if(col == row+2)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+3)
	{
	  return 1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return -3.0/2.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 8)
    {
      if(col ==row)
	{
	  return 7.0/2.0*x;
	}
      else if(col == row-1)
	{
	  return -9.0/2.0*x;
	}
      else if(col == row-2)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row-3)
	{
	  return -1.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return 3.0/2.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }
  else
    {
      if(col == row-1)
	{
	  return -x;
	}
      else if(col == row-2)
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+1)
	{
	  return x;
	}
      else if(col == row+2)
	{
	  return 1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return 3.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return -3.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return 3.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return -3.0/4.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
}

double FinDifGrid::g2value(int row, int col)
{
  double x = (x_dim * x_dim * x_dim)/(x_frct*x_frct*x_frct);

  if(gsl_vector_int_get(typevector,row) == 1)
    {
      if(col == row)
	{
	  return -x;
	}
      else if(col == row+1)
	{
	  return 3.0*x;
	}
      else if(col == row+2)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return -3.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+2)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,4))
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;

	}
    }
  else if(gsl_vector_int_get(typevector,row) == 2)
    {
      if(col == row)
	{
	  return -x;
	}
      else if(col == row-1)
	{
	  return 3.0*x;
	}
      else if(col == row-2)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return -3.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-2)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,4))
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;

	}
    }
  else if(gsl_vector_int_get(typevector,row) == 3)
    {
      if(col == row)
	{
	  return x;
	}
      else if(col == row+1)
	{
	  return -3.0*x;
	}
      else if(col == row+2)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return 3.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+2)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,9))
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;

	}
    }
  else if(gsl_vector_int_get(typevector,row) == 4)
    {
      if(col == row)
	{
	  return x;
	}
      else if(col == row-1)
	{
	  return -3.0*x;
	}
      else if(col == row-2)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return 3.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-2)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,9))
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;

	}
    }


  else if(gsl_vector_int_get(typevector,row) == 5)
    {
      if(col == row)
	{
	  return 7.0/2.0*x;
	}
      else if(col == row-1)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row +1)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -9.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,4))
	{
	  return -1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 6)
    {
      if(col == row)
	{
	  return -7.0/2.0*x;
	}
      else if(col == row-1)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row +1)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return 9.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,9))
	{
	  return 1.0/2.0*x;
	}
      else
	{
	  return 0.0;
	}
    }

  else if(gsl_vector_int_get(typevector,row) == 7)
    {

      if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return 5.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+2)
	{
	  return 3.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -5.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+2)
	{
	  return -3.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return 1.0/4.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }
  else if(gsl_vector_int_get(typevector,row) == 8)
    {

      if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return 5.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return -3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-2)
	{
	  return 3.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return -5.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return 3.0/2.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-2)
	{
	  return -3.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return 1.0/4.0*x;
	}
      else 
	{
	  return 0.0;
	}
    }

  else
    {
      if(col == row+gsl_matrix_int_get(distancetable,row,6))
	{
	  return -x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)-1)
	{
	  return 3.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,6)+1)
	{
	  return 3.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,5))
	{
	  return -1.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7))
	{
	  return x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)-1)
	{
	  return -3.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,7)+1)
	{
	  return -3.0/4.0*x;
	}
      else if(col == row+gsl_matrix_int_get(distancetable,row,8))
	{
	  return 1.0/4.0*x;
	}
      else
	{
	  return 0.0;
	}
    }
}

void FinDifGrid::writea1(gsl_matrix *output)
{
  for(int i = 0; i < fieldpixels; i++)
    {
      for(int j = 0; j < fieldpixels; j++)
	{
	  gsl_matrix_set(output,i,j,a1value(i,j));
	}
    }
}
void FinDifGrid::writea2(gsl_matrix *output)
{
  for(int i = 0; i < fieldpixels; i++)
    {
      for(int j = 0; j < fieldpixels; j++)
	{
	  gsl_matrix_set(output,i,j,a2value(i,j));
	}
    }
}
void FinDifGrid::writes1(gsl_matrix *output)
{
  for(int i = 0; i < fieldpixels; i++)
    {
      for(int j = 0; j < fieldpixels; j++)
	{
	  gsl_matrix_set(output,i,j,s1value(i,j));
	}
    }
}
void FinDifGrid::writes2(gsl_matrix *output)
{
  for(int i = 0; i < fieldpixels; i++)
    {
      for(int j = 0; j < fieldpixels; j++)
	{
	  gsl_matrix_set(output,i,j,s2value(i,j));
	}
    }
}
void FinDifGrid::writec(gsl_matrix *output)
{
  for(int i = 0; i < fieldpixels; i++)
    {
      for(int j = 0; j < fieldpixels; j++)
	{
	  gsl_matrix_set(output,i,j,cvalue(i,j));
	}
    }
}
void FinDifGrid::writef1(gsl_matrix *output)
{
  for(int i = 0; i < fieldpixels; i++)
    {
      for(int j = 0; j < fieldpixels; j++)
	{
	  gsl_matrix_set(output,i,j,f1value(i,j));
	}
    }
}
void FinDifGrid::writef2(gsl_matrix *output)
{
  for(int i = 0; i < fieldpixels; i++)
    {
      for(int j = 0; j < fieldpixels; j++)
	{
	  gsl_matrix_set(output,i,j,f2value(i,j));
	}
    }
}
void FinDifGrid::writeg1(gsl_matrix *output)
{
  for(int i = 0; i < fieldpixels; i++)
    {
      for(int j = 0; j < fieldpixels; j++)
	{
	  gsl_matrix_set(output,i,j,g1value(i,j));
	}
    }
}
void FinDifGrid::writeg2(gsl_matrix *output)
{
  for(int i = 0; i < fieldpixels; i++)
    {
      for(int j = 0; j < fieldpixels; j++)
	{
	  gsl_matrix_set(output,i,j,g2value(i,j));
	}
    }
}

int FinDifGrid::firstorderscheme(int row, int pos)
{
  if(pos == 0 && gsl_matrix_int_get(distancetable,row,5) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,5);
    }
  else if(pos == 1 && gsl_matrix_int_get(distancetable,row,6) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,6);
    }
  else if(pos == 2)
    {
      return row-2;
    }
  else if(pos == 3)
    {
      return row-1;
    }
  else if (pos == 4)
    {
      return row;
    }
  else if(pos == 5)
    {
      return row+1;
    }
  else if(pos == 6)
    {
      return row+2;
    }
  else if(pos == 7 && gsl_matrix_int_get(distancetable,row,7) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,7);
    }
  else if(pos == 8 && gsl_matrix_int_get(distancetable,row,8) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,8);
    }
  else
    {
      return -1;
    }
}

int FinDifGrid::secondorderscheme(int row, int pos)
{
  if(pos == 0 && gsl_matrix_int_get(distancetable,row,5) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,5);
    }
  else if(pos == 1 && gsl_matrix_int_get(distancetable,row,6) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,6)-1;
    }
  else if(pos == 2 && gsl_matrix_int_get(distancetable,row,6) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,6);
    }
  else if(pos == 3 && gsl_matrix_int_get(distancetable,row,6) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,6)+1;
    }
  else if(pos == 4)
    {
      return row-2;
    }
  else if(pos == 5)
    {
      return row-1;
    }
  else if (pos == 6)
    {
      return row;
    }
  else if(pos == 7)
    {
      return row+1;
    }
  else if(pos == 8)
    {
      return row+2;
    }
  else if(pos == 9 && gsl_matrix_int_get(distancetable,row,7) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,7)-1;
    }
  else if(pos == 10 && gsl_matrix_int_get(distancetable,row,7) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,7);
    }
  else if(pos == 11 && gsl_matrix_int_get(distancetable,row,7) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,7)+1;
    }
  else if(pos == 12 && gsl_matrix_int_get(distancetable,row,8) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,8);
    }
  else
    {
      return -1;
    }
}

int FinDifGrid::thirdorderscheme(int row, int pos)
{
  if(pos == 0 && gsl_matrix_int_get(distancetable,row,4) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,4);
    }
  else if(pos == 1 && gsl_matrix_int_get(distancetable,row,5) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,5)-1;
    }
  else if(pos == 2 && gsl_matrix_int_get(distancetable,row,5) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,5);
    }
  else if(pos == 3 && gsl_matrix_int_get(distancetable,row,5) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,5)+1;
    }
  else if(pos == 4 && gsl_matrix_int_get(distancetable,row,6) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,6)-2;
    }
  else if(pos == 5 && gsl_matrix_int_get(distancetable,row,6) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,6)-1;
    }
  else if(pos == 6 && gsl_matrix_int_get(distancetable,row,6) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,6);
    }
  else if(pos == 7 && gsl_matrix_int_get(distancetable,row,6) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,6)+1;
    }
  else if(pos == 8 && gsl_matrix_int_get(distancetable,row,6) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,6)+2;
    }
  else if(pos == 9)
    {
      return row-3;
    }
  else if(pos == 10)
    {
      return row-2;
    }
  else if(pos == 11)
    {
      return row-1;
    }
  else if (pos == 12)
    {
      return row;
    }
  else if(pos == 13)
    {
      return row+1;
    }
  else if(pos == 14)
    {
      return row+2;
    }
  else if(pos == 15)
    {
      return row+3;
    }
  else if(pos == 16 && gsl_matrix_int_get(distancetable,row,7) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,7)-2;
    }
  else if(pos == 17 && gsl_matrix_int_get(distancetable,row,7) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,7)-1;
    }
  else if(pos == 18 && gsl_matrix_int_get(distancetable,row,7) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,7);
    }
  else if(pos == 19 && gsl_matrix_int_get(distancetable,row,7) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,7)+1;
    }
  else if(pos == 20 && gsl_matrix_int_get(distancetable,row,7) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,7)+2;
    }
  else if(pos == 21 && gsl_matrix_int_get(distancetable,row,8) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,8)-1;
    }
  else if(pos == 22 && gsl_matrix_int_get(distancetable,row,8) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,8);
    }
  else if(pos == 23 && gsl_matrix_int_get(distancetable,row,8) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,8)+1;
    }
  else if(pos == 24 && gsl_matrix_int_get(distancetable,row,9) != 0)
    {
      return row + gsl_matrix_int_get(distancetable,row,9);
    }
  else
    {
      return -1;
    }
}
  




void FinDifGrid::a1multvec(gsl_vector *input, gsl_vector *output)
{
  double result;

  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < fieldpixels; j++)
	{
	  result += a1value(i,j)*gsl_vector_get(input,j);
	}
      gsl_vector_set(output,i,result);
    }
}

void FinDifGrid::a1multvec_fast(gsl_vector *input, gsl_vector *output)
{
  double result;
  int cleverindex;
  
  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < 9; j++)
	{
	  cleverindex = firstorderscheme(i,j);
	  if(cleverindex >= 0 && cleverindex < fieldpixels)
	    {
	    result += a1value(i,cleverindex)*gsl_vector_get(input,cleverindex);
	    }
	}
      gsl_vector_set(output,i,result);
    }
}
void FinDifGrid::a2multvec(gsl_vector *input, gsl_vector *output)
{
  double result;

  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < fieldpixels; j++)
	{
	    result += a2value(i,j)*gsl_vector_get(input,j);
	}
      gsl_vector_set(output,i,result);
    }
}

void FinDifGrid::a2multvec_fast(gsl_vector *input, gsl_vector *output)
{
  double result;
  int cleverindex;
  
  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < 9; j++)
	{
	  cleverindex = firstorderscheme(i,j);
	  if(cleverindex >= 0 && cleverindex < fieldpixels)
	    {
	      result += a2value(i,cleverindex)*gsl_vector_get(input,cleverindex);
	    }
	}
      gsl_vector_set(output,i,result);
    }
}

void FinDifGrid::s1multvec(gsl_vector *input, gsl_vector *output)
{
  double result;

  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < fieldpixels; j++)
	{
	  result += s1value(i,j)*gsl_vector_get(input,j);
	}
      gsl_vector_set(output,i,result);
    }
}

void FinDifGrid::s1multvec_fast(gsl_vector *input, gsl_vector *output)
{
  double result;
  int cleverindex;

  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < 13; j++)
	{
	  cleverindex = secondorderscheme(i,j);
	  if(cleverindex >= 0 && cleverindex < fieldpixels)
	    {
	      result += s1value(i,cleverindex)*gsl_vector_get(input,cleverindex);
	    }
	}
      gsl_vector_set(output,i,result);
    }
}
void FinDifGrid::s2multvec(gsl_vector *input, gsl_vector *output)
{
  double result;

  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < fieldpixels; j++)
	{
	  result += s2value(i,j)*gsl_vector_get(input,j);
	}
      gsl_vector_set(output,i,result);
    }
}

void FinDifGrid::s2multvec_fast(gsl_vector *input, gsl_vector *output)
{
  double result;
  int cleverindex;
  
    for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < 13; j++)
	{
	  cleverindex = secondorderscheme(i,j);
	  if(cleverindex >= 0 && cleverindex < fieldpixels)
	    {
	      result += s2value(i,cleverindex)*gsl_vector_get(input,cleverindex);
	    }
	}
      gsl_vector_set(output,i,result);
    }
}

void FinDifGrid::cmultvec(gsl_vector *input, gsl_vector *output)
{
  double result;

  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < fieldpixels; j++)
	{
	  result += cvalue(i,j)*gsl_vector_get(input,j);
	}
      gsl_vector_set(output,i,result);
    }
}

void FinDifGrid::cmultvec_fast(gsl_vector *input, gsl_vector *output)
{
  double result;
  int cleverindex;
  
  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < 13; j++)
	{
	  cleverindex = secondorderscheme(i,j);
	  if(cleverindex >= 0 && cleverindex < fieldpixels)
	    {
	      result += cvalue(i,cleverindex)*gsl_vector_get(input,cleverindex);
	    }
	}
      gsl_vector_set(output,i,result);
    }
}

void FinDifGrid::f1multvec(gsl_vector *input, gsl_vector *output)
{
  double result;

  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < fieldpixels; j++)
	{
	  result += f1value(i,j)*gsl_vector_get(input,j);
	}
      gsl_vector_set(output,i,result);
    }
}

void FinDifGrid::f1multvec_fast(gsl_vector *input, gsl_vector *output)
{
  double result;
  int cleverindex;
  
  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < 25; j++)
	{
	  cleverindex = thirdorderscheme(i,j);
	  if(cleverindex >= 0 && cleverindex < fieldpixels)
	    {
	      result += f1value(i,cleverindex)*gsl_vector_get(input,cleverindex);
	    }
	}
      gsl_vector_set(output,i,result);
    }
}
void FinDifGrid::f2multvec(gsl_vector *input, gsl_vector *output)
{
  double result;

  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < fieldpixels; j++)
	{
	  result += f2value(i,j)*gsl_vector_get(input,j);
	}
      gsl_vector_set(output,i,result);
    }
}

void FinDifGrid::f2multvec_fast(gsl_vector *input, gsl_vector *output)
{
  double result;
  int cleverindex;

  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < 25; j++)
	{
	  cleverindex = thirdorderscheme(i,j);
	  if(cleverindex >= 0 && cleverindex < fieldpixels)
	    {
	      result += f2value(i,cleverindex)*gsl_vector_get(input,cleverindex);
	    }
	}
      gsl_vector_set(output,i,result);
    }
}
void FinDifGrid::g1multvec(gsl_vector *input, gsl_vector *output)
{
  double result;

  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < fieldpixels; j++)
	{
	  result += g1value(i,j)*gsl_vector_get(input,j);
	}
      gsl_vector_set(output,i,result);
    }
}

void FinDifGrid::g1multvec_fast(gsl_vector *input, gsl_vector *output)
{
  double result;
  int cleverindex;
  
  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < 25; j++)
	{
	  cleverindex = thirdorderscheme(i,j);
	  if(cleverindex >= 0 && cleverindex < fieldpixels)
	    {
	      result += g1value(i,cleverindex)*gsl_vector_get(input,cleverindex);
	    }
	}
      gsl_vector_set(output,i,result);
    }
}
void FinDifGrid::g2multvec(gsl_vector *input, gsl_vector *output)
{
  double result;

  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < fieldpixels; j++)
	{
	  result += g2value(i,j)*gsl_vector_get(input,j);
	}
      gsl_vector_set(output,i,result);
    }
}

void FinDifGrid::g2multvec_fast(gsl_vector *input, gsl_vector *output)
{
  double result;
  int cleverindex;
  
  for(int i = 0; i < fieldpixels; i++)
    {
      result = 0.0;
      for(int j = 0; j < 25; j++)
	{
	  cleverindex = thirdorderscheme(i,j);
	  if(cleverindex >= 0 && cleverindex < fieldpixels)
	    {
	      result += g2value(i,cleverindex)*gsl_vector_get(input,cleverindex);
	    }
	}
      gsl_vector_set(output,i,result);
    }
}

double FinDifGrid::ShearBlk(int l, int k, gsl_vector *ellip1, gsl_vector *ellip2, gsl_matrix *F1, gsl_matrix *F2,string selection)
{

  double result1 = 0.0;
  double result2 = 0.0;
  int kis;
  int kjs;
  int lis;
  int ljs;

  if(selection == "reduced shear")
    {

      for(int i = 0; i < 13; i++)
	{
	  for(int j = 0; j < 13; j++)
	    {
	      kis = secondorderscheme(k,i);
	      ljs = secondorderscheme(l,j);
	      
	      if(kis >= 0 && kis < fieldpixels && ljs >= 0 && ljs < fieldpixels)
		{
		  result1 += gsl_matrix_get(F1,kis,ljs)*(gsl_vector_get(ellip1,kis)*gsl_vector_get(ellip1,ljs)*cvalue(kis,k)*cvalue(ljs,l) + gsl_vector_get(ellip1,kis)*cvalue(kis,k)*s1value(ljs,l)+ s1value(kis,k)*s1value(ljs,l))+ gsl_matrix_get(F2,kis,ljs)*(gsl_vector_get(ellip2,kis)*gsl_vector_get(ellip2,ljs)*cvalue(kis,k)*cvalue(ljs,l)+ gsl_vector_get(ellip2,kis)*cvalue(kis,k)*s2value(ljs,l)+ s2value(kis,k)*s2value(ljs,l));
		}
	    }
	}
      for(int i = 0; i < 13; i++)
	{
	  for(int j = 0; j < 13; j++)
	    {
	      kjs = secondorderscheme(k,j);
	      lis = secondorderscheme(l,i);
	      
	      if(kjs >= 0 && kjs < fieldpixels && lis >= 0 && lis < fieldpixels)
		{
		  result2 += gsl_matrix_get(F1,lis,kjs)*gsl_vector_get(ellip1,lis)*s1value(kjs,k)*cvalue(lis,l) + gsl_matrix_get(F2,lis,kjs)*gsl_vector_get(ellip2,lis)*s2value(kjs,k)*cvalue(lis,l);
		}
	    }
	}
      return 2.0*(result1+result2);
    }

  else if(selection == "shear")
    {
      for(int i = 0; i < 13; i++)
	{
	  for(int j = 0; j < 13; j++)
	    {
	      kis = secondorderscheme(k,i);
	      ljs = secondorderscheme(l,j);
	      
	      if(kis >= 0 && kis < fieldpixels && ljs >= 0 && ljs < fieldpixels)
		{
		  result1 += gsl_matrix_get(F1,kis,ljs)*s1value(kis,k)*s1value(ljs,l)+gsl_matrix_get(F2,kis,ljs)*s2value(kis,k)*s2value(ljs,l);
		}
	    }
	}
      return 2.0*result1;
    }
  else
    {
      throw invalid_argument("Selection not valid for ShearBlk");
    }
}

double FinDifGrid::FlexionBlk(int l, int k,gsl_matrix *F1, gsl_matrix *F2,gsl_matrix *G1, gsl_matrix *G2)
{

  double result1 = 0.0;
  int kis;
  int ljs;

  for(int i = 0; i < 24; i++)
    {
      for(int j = 0; j < 24; j++)
	{
	  kis = thirdorderscheme(k,i);
	  ljs = thirdorderscheme(l,j);
	  
	  if(kis >= 0 && kis < fieldpixels && ljs >= 0 && ljs < fieldpixels)
	    {
	      result1 += gsl_matrix_get(F1,kis,ljs)*f1value(kis,k)*f1value(ljs,l)+gsl_matrix_get(F2,kis,ljs)*f2value(kis,k)*f2value(ljs,l)+gsl_matrix_get(G1,kis,ljs)*g1value(kis,k)*g1value(ljs,l)+gsl_matrix_get(G2,kis,ljs)*g2value(kis,k)*g2value(ljs,l);
	    }
	}
    }
  return 2.0*result1;
  
}

double FinDifGrid::CcurveBlk(int l, int k, gsl_vector_int *ccurve, gsl_vector *factors)
{

  double result = 0.0;
  int kis;

  for(int i = 0; i < 13; i++)
    {
      kis = secondorderscheme(k,i);
      if(kis >= 0 && kis < fieldpixels)
	{
	  if(gsl_vector_int_get(ccurve,kis) != 0)
	    {
	      result += gsl_vector_get(factors,kis)*(cvalue(kis,k)*cvalue(kis,l) - s1value(kis,k)*s1value(kis,l) - s2value(kis,k)*s2value(kis,l));
	    }
	}
    }
  return result;
}

double FinDifGrid::MsystemBlk(int l, int k, gsl_matrix *msysteminfo)
{

  double meanpos1;
  double meanpos2;
  double interresult1;
  double interresult2;
  double nofuture1 = 0.0;
  double nofuture2 = 0.0;


  for(int ms = 0; ms < msysteminfo->size2; ms++)
    {
      meanpos1 = 0.0;
      meanpos2 = 0.0;
      for(int n = 0; n < (int) gsl_matrix_get(msysteminfo,0,ms); n++)
	{
	  for(int m = 0; m < (int) gsl_matrix_get(msysteminfo,0,ms); m++)
	    {
	      meanpos1 += a1value((int) gsl_matrix_get(msysteminfo,12+m,ms)-1,k)*a1value((int) gsl_matrix_get(msysteminfo,12+n,ms)-1,l);
	      meanpos2 += a2value((int) gsl_matrix_get(msysteminfo,12+m,ms)-1,k)*a2value((int) gsl_matrix_get(msysteminfo,12+n,ms)-1,l);
	    }
	}
      meanpos1 = 2.0*meanpos1;
      meanpos2 = 2.0*meanpos2;
      for(int a = 0; a < gsl_matrix_get(msysteminfo,0,ms); a++)
	{
	  interresult1 = 0.0;
	  interresult2 = 0.0;
	  
	  for(int n = 0; n < (int) gsl_matrix_get(msysteminfo,0,ms); n++)
	    {
	      interresult1 += (a1value((int) gsl_matrix_get(msysteminfo,12+n,ms)-1,k)*a1value((int) gsl_matrix_get(msysteminfo,12+n,ms)-1,l) - a1value((int) gsl_matrix_get(msysteminfo,12+a,ms)-1,k)*a1value((int) gsl_matrix_get(msysteminfo,12+n,ms)-1,l));
	      interresult2 += (a2value((int) gsl_matrix_get(msysteminfo,12+n,ms)-1,k)*a2value((int) gsl_matrix_get(msysteminfo,12+n,ms)-1,l) - a2value((int) gsl_matrix_get(msysteminfo,12+a,ms)-1,k)*a2value((int) gsl_matrix_get(msysteminfo,12+n,ms)-1,l));
	    }
	  
	  interresult1 = interresult1/gsl_matrix_get(msysteminfo,0,ms);
	  interresult2 = interresult2/gsl_matrix_get(msysteminfo,0,ms);
	  nofuture1 += a1value((int) gsl_matrix_get(msysteminfo,12+a,ms)-1,k)*a1value((int) gsl_matrix_get(msysteminfo,12+a,ms)-1,l) + interresult1;
	  nofuture2 += a2value((int) gsl_matrix_get(msysteminfo,12+a,ms)-1,k)*a2value((int) gsl_matrix_get(msysteminfo,12+a,ms)-1,l) + interresult2;
	}
      nofuture1 = nofuture1*2.0/(gsl_matrix_get(msysteminfo,17,ms)*gsl_matrix_get(msysteminfo,17,ms));
      nofuture2 = nofuture2*2.0/(gsl_matrix_get(msysteminfo,17,ms)*gsl_matrix_get(msysteminfo,17,ms));
      nofuture1 += meanpos1;
      nofuture2 += meanpos2;
    }
  return nofuture1+nofuture2;
}

double FinDifGrid::RegularisationBlk(int l, int k, double eta1, double eta2, string selection)
{

  double result = 0.0;
  double result2 = 0.0;
  int kis;

  if(selection == "convergence")
    {
      for(int i = 0; i < 13; i++)
	{
	  kis = secondorderscheme(l,i);
	  if(kis >= 0 && kis < fieldpixels)
	    {
	      result += cvalue(kis,l)*cvalue(kis,k);
	    }
	}
      return eta1*result;
    }
  else if(selection == "convshear")
    {
      for(int i = 0; i < 13; i++)
	{
	  kis = secondorderscheme(l,i);
	  if(kis >= 0 && kis < fieldpixels)
	    {
	      result += cvalue(kis,l)*cvalue(kis,k) + s1value(kis,l)*s1value(kis,k) + s2value(kis,l)*s2value(kis,k) ;
	    }
	}
      return eta1*result;
    }
  else if(selection == "flexion")
    {
      for(int i = 0; i < 25; i++)
	{
	  kis = thirdorderscheme(l,i);
	  if(kis >= 0 && kis < fieldpixels)
	    {
	      result += f1value(kis,l)*f1value(kis,k) + f2value(kis,l)*f2value(kis,k) + g1value(kis,l)*g1value(kis,k) + g2value(kis,l)*g2value(kis,k) ;
	    }
	}
      return eta2*result;
    }
  else if(selection == "convshearflex")
    {
      for(int i = 0; i < 13; i++)
	{
	  kis = secondorderscheme(l,i);
	  if(kis >= 0 && kis < fieldpixels)
	    {
	      result += cvalue(kis,l)*cvalue(kis,k) + s1value(kis,l)*s1value(kis,k) + s2value(kis,l)*s2value(kis,k) ;
	    }
	}
      for(int i = 0; i < 25; i++)
	{
	  kis = thirdorderscheme(l,i);
	  if(kis >= 0 && kis < fieldpixels)
	    {
	      result2 += f1value(kis,l)*f1value(kis,k) + f2value(kis,l)*f2value(kis,k) + g1value(kis,l)*g1value(kis,k) + g2value(kis,l)*g2value(kis,k) ;
	    }
	}
      return eta1*result + eta2*result2;
    }
  else
    {
      throw invalid_argument("Selection not valid for RegularisationBlk");
    }
}

double FinDifGrid::ShearVl(int l, gsl_vector *ellip1, gsl_vector * ellip2, gsl_matrix *F1, gsl_matrix * F2, gsl_vector *redshift, string selection)
{

  double result1 = 0.0;
  double result2 = 0.0;
  int lis;
  int ljs;

  if(selection == "reduced shear")
    {

      for(int i = 0; i < fieldpixels; i++)
	{
	  for(int j = 0; j < 13; j++)
	    {
	      ljs = secondorderscheme(l,j);
	      if(ljs >= 0 && ljs < fieldpixels)
		{
		  result1 += gsl_matrix_get(F1,i,ljs)*(gsl_vector_get(ellip1,i)*gsl_vector_get(ellip1,ljs)*gsl_vector_get(redshift,ljs)*cvalue(ljs,l)+gsl_vector_get(ellip1,i)*gsl_vector_get(redshift,ljs)*s1value(ljs,l))+gsl_matrix_get(F2,i,ljs)*(gsl_vector_get(ellip2,i)*gsl_vector_get(ellip2,ljs)*gsl_vector_get(redshift,ljs)*cvalue(ljs,l)+gsl_vector_get(ellip2,i)*gsl_vector_get(redshift,ljs)*s2value(ljs,l));
		}
	    }
	}
      for(int i = 0; i < 13; i++)
	{
	  for(int j = 0; j < fieldpixels; j++)
	    {
	      lis = secondorderscheme(l,i);
	      if(lis >= 0 && lis < fieldpixels)
		{
		  result2 += gsl_matrix_get(F1,lis,j)*(gsl_vector_get(ellip1,lis)*gsl_vector_get(ellip1,j)*gsl_vector_get(redshift,lis)*cvalue(lis,l)+gsl_vector_get(ellip1,j)*gsl_vector_get(redshift,lis)*s1value(lis,l))+gsl_matrix_get(F2,lis,j)*(gsl_vector_get(ellip2,lis)*gsl_vector_get(ellip2,j)*gsl_vector_get(redshift,lis)*cvalue(lis,l)+gsl_vector_get(ellip2,j)*gsl_vector_get(redshift,lis)*s2value(lis,l));
		}
	    }
	}
      return result1+result2;
    }
  else if(selection == "shear")
    {
      for(int i = 0; i < fieldpixels; i++)
	{
	  for(int j = 0; j < 13; j++)
	    {
	      ljs = secondorderscheme(l,j);
	      if(ljs >= 0 && ljs < fieldpixels)
		{
		  result1 += gsl_matrix_get(F1,i,ljs)*gsl_vector_get(ellip1,i)*gsl_vector_get(redshift,ljs)*s1value(ljs,l) + gsl_matrix_get(F2,i,ljs)*gsl_vector_get(ellip2,i)*gsl_vector_get(redshift,ljs)*s2value(ljs,l); 
		}
	    }
	}
      for(int i = 0; i < 13; i++)
	{
	  for(int j = 0; j < fieldpixels; j++)
	    {
	      lis = secondorderscheme(l,i);
	      if(lis >= 0 && lis < fieldpixels)
		{
		  result2 += gsl_matrix_get(F1,lis,j)*gsl_vector_get(ellip1,j)*gsl_vector_get(redshift,lis)*s1value(lis,l) + gsl_matrix_get(F2,lis,j)*gsl_vector_get(ellip2,j)*gsl_vector_get(redshift,lis)*s2value(lis,l);
		}
	    }
	}
      return result1+result2;
    }
  else
    {
      throw invalid_argument("Selection not valid for ShearVl");
    }

}


double FinDifGrid::FlexionVl(int l, gsl_vector *f1, gsl_vector *f2, gsl_vector *g1, gsl_vector *g2, gsl_matrix *F1, gsl_matrix *F2, gsl_matrix *G1, gsl_matrix *G2, gsl_vector *redshift)
{

  double result = 0.0;
  int ljs;

  for(int i = 0; i < fieldpixels; i++)
    {
      for(int j = 0; j < 25; j++)
	{
	  ljs = thirdorderscheme(l,j);
	  if(ljs >= 0 && ljs < fieldpixels)
	    {
	      result += gsl_vector_get(redshift,ljs)*(gsl_matrix_get(F1,i,ljs)*gsl_vector_get(f1,i)*f1value(ljs,l) + gsl_matrix_get(F2,i,ljs)*gsl_vector_get(f2,i)*f2value(ljs,l) + gsl_matrix_get(G1,i,ljs)*gsl_vector_get(g1,i)*g1value(ljs,l) + gsl_matrix_get(G2,i,ljs)*gsl_vector_get(g2,i)*g2value(ljs,l));
	    }
	}
    }
  return 2.0*result;
}

double FinDifGrid::CcurveVl(int l, gsl_vector_int *ccurve,gsl_vector *factors)
{

  double result = 0.0;
  int lis;

  for(int i = 0; i < 13; i++)
    {
      lis = secondorderscheme(l,i);
      if(lis >= 0 && lis < fieldpixels)
	{
	  if(gsl_vector_int_get(ccurve,lis) != 0)
	    {
	      result += gsl_vector_get(factors,lis)*cvalue(lis,l);
	    }
	}
    }
  return result;
}

double FinDifGrid::MsystemVl(int l, gsl_matrix *msysteminfo)
{
  double meanpos1;
  double meanpos2;
  double interresult1;
  double interresult2;
  double nofuture1 = 0.0;
  double nofuture2 = 0.0;

  for(int ms = 0; ms < msysteminfo->size2; ms++)
    {
      meanpos1 = 0.0;
      meanpos2 = 0.0;
      
      for(int n = 0; n < (int) gsl_matrix_get(msysteminfo,0,ms); n++)
	{
	  for(int m = 0; m < (int) gsl_matrix_get(msysteminfo,0,ms); m++)
	    {
	      meanpos1 += gsl_matrix_get(msysteminfo,2+2*m,ms)*a1value((int) gsl_matrix_get(msysteminfo,12+n,ms)-1,l);
	      meanpos2 += gsl_matrix_get(msysteminfo,3+2*m,ms)*a2value((int) gsl_matrix_get(msysteminfo,12+n,ms)-1,l);
	    }
	}
      meanpos1 = 2.0*meanpos1;
      meanpos2 = 2.0*meanpos2;
      for(int a = 0; a < (int) gsl_matrix_get(msysteminfo,0,ms); a++)
	{
	  interresult1 = 0.0;
	  interresult2 = 0.0;
	  
	  for(int n = 0; n < (int) gsl_matrix_get(msysteminfo,0,ms); n++)
	    {
	      interresult1 += (gsl_matrix_get(msysteminfo,2+2*a,ms)*a1value((int) gsl_matrix_get(msysteminfo,12+n,ms)-1,l) + gsl_matrix_get(msysteminfo,2+2*n,ms)*a1value((int) gsl_matrix_get(msysteminfo,12+n,ms)-1,l));
	      interresult2 += (gsl_matrix_get(msysteminfo,3+2*a,ms)*a2value((int) gsl_matrix_get(msysteminfo,12+n,ms)-1,l) + gsl_matrix_get(msysteminfo,3+2*n,ms)*a2value((int) gsl_matrix_get(msysteminfo,12+n,ms)-1,l));
	    }
	  interresult1 = interresult1/gsl_matrix_get(msysteminfo,0,ms);
	  interresult2 = interresult2/gsl_matrix_get(msysteminfo,0,ms);

	  nofuture1 += interresult1 -  gsl_matrix_get(msysteminfo,2+2*a,ms)*a1value((int) gsl_matrix_get(msysteminfo,12+a,ms)-1,l);
	  nofuture2 += interresult2 -  gsl_matrix_get(msysteminfo,3+2*a,ms)*a2value((int) gsl_matrix_get(msysteminfo,12+a,ms)-1,l);
	}
      nofuture1 = nofuture1*2.0/(gsl_matrix_get(msysteminfo,17,ms)*gsl_matrix_get(msysteminfo,17,ms));
      nofuture2 = nofuture2*2.0/(gsl_matrix_get(msysteminfo,17,ms)*gsl_matrix_get(msysteminfo,17,ms));
	  nofuture1 += meanpos1;
	  nofuture2 += meanpos2;
    }

  return nofuture1+nofuture2;
}

double FinDifGrid::RegularisationVl(int l, gsl_vector *conv, gsl_vector *shear1, gsl_vector *shear2, gsl_vector *f1, gsl_vector *f2, gsl_vector *g1, gsl_vector *g2, double eta1, double eta2, string selection)
{

  double result = 0.0;
  double result2 = 0.0;
  int lis;

  if(selection == "convergence")
    {
      for(int i = 0; i < 13; i++)
	{
	  lis = secondorderscheme(l,i);
	  if(lis >= 0 && lis < fieldpixels)
	    {
	      result += gsl_vector_get(conv,lis)*cvalue(lis,l);
	    }
	}
      return eta1*result;
    }
  else if(selection == "convshear")
    {
      for(int i = 0; i < 13; i++)
	{
	  lis = secondorderscheme(l,i);
	  if(lis >= 0 && lis < fieldpixels)
	    {
	      result += gsl_vector_get(conv,lis)*cvalue(lis,l) + gsl_vector_get(shear1,lis)*s1value(lis,l) + gsl_vector_get(shear2,lis)*s2value(lis,l);
	    }
	}
      return eta1*result;
    }
  else if(selection == "flexion")
    {
      for(int i = 0; i < 25; i++)
	{
	  lis = thirdorderscheme(l,i);
	  if(lis >= 0 && lis < fieldpixels)
	    {
	      result += gsl_vector_get(f1,lis)*f1value(lis,l) + gsl_vector_get(f2,lis)*f2value(lis,l) + gsl_vector_get(g1,lis)*g1value(lis,l) + gsl_vector_get(g2,lis)*g2value(lis,l);
	    }
	}
      return eta2*result;
    }
  else if(selection == "convshearflex")
    {
      for(int i = 0; i < 13; i++)
	{
	  lis = secondorderscheme(l,i);
	  if(lis >= 0 && lis < fieldpixels)
	    {
	      result += gsl_vector_get(conv,lis)*cvalue(lis,l) + gsl_vector_get(shear1,lis)*s1value(lis,l) + gsl_vector_get(shear2,lis)*s2value(lis,l);
	    }
	}
      for(int i = 0; i < 25; i++)
	{
	  lis = thirdorderscheme(l,i);
	  if(lis >= 0 && lis < fieldpixels)
	    {
	      result += gsl_vector_get(f1,lis)*f1value(lis,l) + gsl_vector_get(f2,lis)*f2value(lis,l) + gsl_vector_get(g1,lis)*g1value(lis,l) + gsl_vector_get(g2,lis)*g2value(lis,l);
	    }
	}

      return eta1*result + eta2*result2;
    }
  else
    {
      throw invalid_argument("Selection not valid for RegularisationVl");
    }
}
 





