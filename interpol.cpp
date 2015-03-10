#include "interpol.h"


void cspline_mapinterpol(gsl_matrix *in, gsl_matrix *out)
{
  int xin = in->size2;
  int yin = in->size1;
  int xout = out->size2;
  int yout = out->size1;


  gsl_matrix *intermediate = gsl_matrix_calloc(xout,yin);
  double value;
  gsl_vector *support = gsl_vector_calloc(xin);
  double instep = 1.0/xin;
  double outstep = 1.0/xout;


  //creating the rows of the interpolation

  for(int i = 0; i < xin; i++)
    {
      gsl_vector_set(support,i,instep/2.0+i*instep);
    }



  for(int i = 0; i < yin; i++)
    {
      gsl_interp_accel *acc = gsl_interp_accel_alloc ();
      gsl_spline *spline  = gsl_spline_alloc (gsl_interp_cspline, xin);
      gsl_spline_init (spline, gsl_vector_ptr(support,0),gsl_matrix_ptr(in,i,0),xin);

      for(int j = 0; j < xout; j++)
	{
	  gsl_matrix_set(intermediate,j,i,gsl_spline_eval(spline,outstep/2.0+j*outstep,acc));
	}

      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);
    }

  for(int i= 0; i < xout; i++)
    {
      gsl_interp_accel *acc = gsl_interp_accel_alloc ();
      gsl_spline *spline  = gsl_spline_alloc (gsl_interp_cspline, yin);
      gsl_spline_init (spline, gsl_vector_ptr(support,0),gsl_matrix_ptr(intermediate,i,0),yin);


      for(int j = 0; j < yout; j++)
	{
	  gsl_matrix_set(out,j,i,gsl_spline_eval(spline,outstep/2.0+j*outstep,acc));
	  
	}

      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);
    }

  gsl_matrix_free(intermediate);

}

void cspline_mapinterpol_smooth(gsl_matrix *in, gsl_matrix *out)
{
  int xin = in->size2;
  int yin = in->size1;
  int xout = out->size2;
  int yout = out->size1;

  gsl_matrix *inter;
  gsl_matrix *nofuture = gsl_matrix_calloc(yin,xin);
  int length = xout-xin;
  int yvalue1 = 0;
  int yvalue2 = 0;
  double ratio = (double) xout/yout;
  cout <<"test" <<endl;
  cout <<yin <<endl;
  cout <<xin <<endl;
  gsl_matrix_memcpy(nofuture,in);

  for(int i = 1; i <= length; i++)
    {
      yvalue1 = (int) ceil((xin+i-1)/ratio);
      yvalue2 = (int) ceil((xin+i)/ratio);
      inter = gsl_matrix_calloc(yvalue2,xin+i);
      cspline_mapinterpol(nofuture,inter);
      gsl_matrix_free(nofuture);
      nofuture = gsl_matrix_calloc(yvalue2,xin+i);
      gsl_matrix_memcpy(nofuture,inter);
      gsl_matrix_free(inter);
    }

  cout <<nofuture->size1 <<endl;
  cout <<nofuture->size2 <<endl;
  gsl_matrix_memcpy(out,nofuture);
  gsl_matrix_free(nofuture);
}


void cspline_maskedinterpol_row(gsl_matrix *in, gsl_matrix_int *inmap, gsl_matrix *out)
{

  int xin = in->size2;
  int yin = in->size1;
  int xout = out->size2;
  int yout = out->size1;

  gsl_matrix *intermediate = gsl_matrix_calloc(yin,xout);
  double value;
  double instep = 1.0/xin;
  double outstep = 1.0/xout;
  int counter;
  int index;
  gsl_vector *support;
  gsl_vector *data;
  gsl_vector_int *badrow = gsl_vector_int_calloc(yin);
  gsl_interp_accel *acc;
  gsl_spline *spline;

  for(int i = 0; i < yin; i++)
    {
      counter = 0;
      for(int j = 0; j < xin; j++)
	{
	  if(gsl_matrix_int_get(inmap,i,j) != 1)
	    {
	      counter++;
	    }
	}
      if(counter != 0)
	{
	  data = gsl_vector_calloc(counter);
	  support = gsl_vector_calloc(counter);
	  index = 0;
	  
	  for(int j = 0; j < xin; j++)
	    {
	      if(gsl_matrix_int_get(inmap,i,j) != 1)
		{
		  gsl_vector_set(data,index,gsl_matrix_get(in,i,j));
		  gsl_vector_set(support,index,instep/2.0+j*instep);
		  index++;
		}
	    }
	  
	  acc = gsl_interp_accel_alloc ();
	  spline  = gsl_spline_alloc (gsl_interp_cspline, counter);
	  gsl_spline_init (spline, gsl_vector_ptr(support,0),gsl_vector_ptr(data,0),counter);
	  
	  for(int j = 0; j < xout; j++)
	    {
	      gsl_matrix_set(intermediate,i,j,gsl_spline_eval(spline,outstep/2.0+j*outstep,acc));
	    }
	  gsl_vector_int_set(badrow,i,0);
	  gsl_vector_free(data);
	  gsl_vector_free(support);
	  gsl_spline_free(spline);
	  gsl_interp_accel_free(acc);
	}

      else
	{
	  gsl_vector_int_set(badrow,i,1);
	}

    }






  for(int i= 0; i < xout; i++)
    {
      counter = 0;
      for(int j = 0; j < yin; j++)
	{
	  if(gsl_vector_int_get(badrow,j) != 1)
	    {
	      counter++;
	    }
	}
      support = gsl_vector_calloc(counter);
      data = gsl_vector_calloc(counter);

      index = 0;
      for(int j = 0; j < yin; j++)
	{
	  if(gsl_vector_int_get(badrow,j) != 1)
	    {
	      gsl_vector_set(data,index,gsl_matrix_get(intermediate,j,i));
	      gsl_vector_set(support,index,instep/2.0+j*instep);
	      index++;
	    }
	}

      acc = gsl_interp_accel_alloc();
      spline  = gsl_spline_alloc (gsl_interp_cspline, counter);
      gsl_spline_init (spline, gsl_vector_ptr(support,0),gsl_vector_ptr(data,0),counter);


      for(int j = 0; j < yout; j++)
	{
	  gsl_matrix_set(out,j,i,gsl_spline_eval(spline,outstep/2.0+j*outstep,acc));
	  
	}

      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);
      gsl_vector_free(data);
      gsl_vector_free(support);

    }

  gsl_matrix_free(intermediate);

}


void cspline_maskedinterpol_col(gsl_matrix *in, gsl_matrix_int *inmap, gsl_matrix *out)
{

  int xin = in->size2;
  int yin = in->size1;
  int xout = out->size2;
  int yout = out->size1;

  gsl_matrix *intermediate = gsl_matrix_calloc(yout,xin);
  double value;
  double instep = 1.0/xin;
  double outstep = 1.0/xout;
  int counter;
  int index;
  gsl_vector *support;
  gsl_vector *data;
  gsl_vector_int *badcol = gsl_vector_int_calloc(xin);
  gsl_interp_accel *acc;
  gsl_spline *spline;

  for(int i = 0; i < xin; i++)
    {
      counter = 0;
      for(int j = 0; j < yin; j++)
	{
	  if(gsl_matrix_int_get(inmap,j,i) != 1)
	    {
	      counter++;
	    }
	}
      if(counter != 0)
	{
	  data = gsl_vector_calloc(counter);
	  support = gsl_vector_calloc(counter);
	  index = 0;
	  
	  for(int j = 0; j < yin; j++)
	    {
	      if(gsl_matrix_int_get(inmap,j,i) != 1)
		{
		  gsl_vector_set(data,index,gsl_matrix_get(in,j,i));
		  gsl_vector_set(support,index,instep/2.0+j*instep);
		  index++;
		}
	    }
	  
	  acc = gsl_interp_accel_alloc ();
	  spline  = gsl_spline_alloc (gsl_interp_cspline, counter);
	  gsl_spline_init (spline, gsl_vector_ptr(support,0),gsl_vector_ptr(data,0),counter);
	  
	  for(int j = 0; j < yout; j++)
	    {
	      gsl_matrix_set(intermediate,j,i,gsl_spline_eval(spline,outstep/2.0+j*outstep,acc));
	    }
	  gsl_vector_int_set(badcol,i,0);
	  gsl_vector_free(data);
	  gsl_vector_free(support);
	  gsl_spline_free(spline);
	  gsl_interp_accel_free(acc);
	}

      else
	{
	  gsl_vector_int_set(badcol,i,1);
	}

    }







  for(int i= 0; i < yout; i++)
    {
      counter = 0;
      for(int j = 0; j < xin; j++)
	{
	  if(gsl_vector_int_get(badcol,j) != 1)
	    {
	      counter++;
	    }
	}
      support = gsl_vector_calloc(counter);
      data = gsl_vector_calloc(counter);

      index = 0;
      for(int j = 0; j < xin; j++)
	{
	  if(gsl_vector_int_get(badcol,j) != 1)
	    {
	      gsl_vector_set(data,index,gsl_matrix_get(intermediate,i,j));
	      gsl_vector_set(support,index,instep/2.0+j*instep);
	      index++;
	    }
	}

      acc = gsl_interp_accel_alloc();
      spline  = gsl_spline_alloc (gsl_interp_cspline, counter);
      gsl_spline_init (spline, gsl_vector_ptr(support,0),gsl_vector_ptr(data,0),counter);


      for(int j = 0; j < xout; j++)
	{
	  gsl_matrix_set(out,i,j,gsl_spline_eval(spline,outstep/2.0+j*outstep,acc));
	  
	}

      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);
      gsl_vector_free(data);
      gsl_vector_free(support);

    }

  gsl_matrix_free(intermediate);

}




void cspline_morecleverinterpol(gsl_matrix *in, gsl_matrix_int *inmask, gsl_matrix *out)
{

  int xin = in->size2;
  int yin = in->size1;
  int xout = out->size2;
  int yout = out->size1;



  //preapring the inputmap, by setting each masked pixel to the value
  //of it's nearesat neighbour

  double distanceact;
  double value = 0.0;
  
  for(int i = 0; i < yin; i++)
    {
      for(int j = 0; j < xin; j++)
	{
	  double distancemin = (int) xin + (int) yin;
	  if(gsl_matrix_int_get(inmask,i,j) == 1)
	    {
	      for(int k = 0; k < yin; k++)
		{
		  for(int l = 0; l < xin; l++)
		    {
		      distanceact = sqrt(pow(i-k,2.0)+pow(j-l,2.0));
		      if(gsl_matrix_int_get(inmask,k,l) != 1 && distanceact < distancemin)
			{
			  distancemin = distanceact;
			  value = gsl_matrix_get(in,k,l);
			}
		    }
		}
	      gsl_matrix_set(in,i,j,value);
	    }
	}
    }
  



  gsl_matrix_int *dummy1 = gsl_matrix_int_calloc(yin,xin);
  gsl_matrix_int *dummy2 = gsl_matrix_int_calloc(yout,xout);
  cspline_cleverinterpol(in,dummy1,out,dummy2);

}




 

void cspline_cleverinterpol(gsl_matrix *in, gsl_matrix_int *inmask, gsl_matrix* out, gsl_matrix_int *outmask)
{
  int xin = in->size2;
  int yin = in->size1;
  int xout = out->size2;
  int yout = out->size1;

  gsl_matrix *samplerow = gsl_matrix_calloc(yout,xout);
  gsl_matrix *samplecol = gsl_matrix_calloc(yout,xout);
  gsl_vector *sample1;
  gsl_vector *rowconv;
  gsl_vector *colconv;
  int counter = 0;
  int index;
  double distance;
  double wmean1 = 0;
  double wmean2 = 0;
  int row = 0;
  int col = 0;

  for(int i = 0; i < yout; i++)
    {
      for(int j = 0; j < xout; j++)
	{
	  if(gsl_matrix_int_get(outmask,i,j) != 1)
	    {
	      counter++;
	    }
	}
    }

  rowconv = gsl_vector_calloc(counter);
  colconv = gsl_vector_calloc(counter);
  sample1 = gsl_vector_calloc(2*counter);

  cspline_maskedinterpol_row(in, inmask, samplerow);
  cspline_maskedinterpol_col(in, inmask, samplecol);


  
  index = 0;
  for(int i = 0; i < yout; i++)
    {
      for(int j = 0; j < xout; j++)
	{
	  if(gsl_matrix_int_get(outmask,i,j) != 1)
	  {
	      gsl_vector_set(sample1,index,gsl_matrix_get(samplerow,i,j));
	      gsl_vector_set(sample1,index+counter,gsl_matrix_get(samplecol,i,j));
	      index++;
	  }
	}
    }
  

  wmean1 = gsl_stats_mean(gsl_vector_ptr(sample1,0),1,2*counter);

  gsl_matrix *map1 = gsl_matrix_calloc(yout,xout);
  gsl_matrix *map2 = gsl_matrix_calloc(yout,xout);
  gsl_matrix_int *map3 = gsl_matrix_int_calloc(yout,xout);

  for(int i = 0; i < yout; i++)
    {
      for(int j = 0; j < xout; j++)
	{

	  gsl_matrix_set(map1,i,j,abs(gsl_matrix_get(samplerow,i,j)-wmean1));
	  gsl_matrix_set(map2,i,j,abs(gsl_matrix_get(samplecol,i,j)-wmean1));

	  if(abs(gsl_matrix_get(samplerow,i,j)-wmean1) <= abs(gsl_matrix_get(samplecol,i,j)-wmean1))
	    {
	      gsl_matrix_set(out,i,j,gsl_matrix_get(samplerow,i,j));
	      row++;
	      gsl_matrix_int_set(map3,i,j,10);
	    }
	  else
	    {
	      gsl_matrix_set(out,i,j,gsl_matrix_get(samplecol,i,j));
	      col++;
	      gsl_matrix_int_set(map3,i,j,20);
	    }
	}
    
    }




}

void cspline_smooth_cleverinterpol(gsl_matrix *in, gsl_matrix_int *inmask, gsl_matrix* out,gsl_matrix_int *outmask)
{
  int xin = in->size2;
  int yin = in->size1;
  int xout = out->size2;
  int yout = out->size1;

  double ratio = xout/yout;
  int yvalue = (int) outmask->size1;
  gsl_matrix *nofuture = gsl_matrix_calloc(yvalue,xin+1);
  cspline_cleverinterpol(in,inmask,nofuture,outmask);
  cspline_mapinterpol_smooth(nofuture,out);
}

void cspline_smooth_morecleverinterpol(gsl_matrix *in, gsl_matrix_int *inmask, gsl_matrix* out)
{
  int xin = in->size2;
  int yin = in->size1;


  int ratio = (int) floor((xin+1)*yin/xin);
  gsl_matrix *nofuture = gsl_matrix_calloc(ratio,xin+1);
  cspline_morecleverinterpol(in,inmask,nofuture);
  cspline_mapinterpol_smooth(nofuture,out);
}













			     












