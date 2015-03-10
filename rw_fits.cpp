#include "rw_fits.h"

using namespace CCfits;

void write_pimg(string filename, gsl_matrix *data)
{

  long naxis = 2;
  int x_dim = data->size2;
  int y_dim = data->size1;
  long naxes[2] = {x_dim, y_dim};
  int dim = x_dim * y_dim;
  std::auto_ptr<FITS> img(0);
  filename = "!"+filename;
  img.reset(new FITS(filename,DOUBLE_IMG,naxis,naxes));

  std::valarray<double> array(dim); 
  for (int i = 0; i < y_dim; i++)
    {
      for(int j = 0; j < x_dim; j++)
	{
	 array[i*x_dim + j] = gsl_matrix_get(data, i, j);
	}
    }
  

  img->pHDU().write(1, dim, array);
}

void write_pimgint(string filename, gsl_matrix_int *data)
{
  //same as above with integer values

  int x_dim = data->size2;
  int y_dim = data->size1;
  long naxis = 2;
  long naxes[2] = {x_dim, y_dim};
  int dim = x_dim * y_dim;
  filename = "!"+filename;

  std::auto_ptr<FITS> img(0);
  img.reset(new FITS(filename,32,naxis,naxes));

  std::valarray<int> array(dim); 
  for (int i = 0; i < y_dim; i++)
    {
      for(int j = 0; j < x_dim; j++)
	{
	  array[i*x_dim + j] = gsl_matrix_int_get(data, i, j);
	}
    }
  
  img->pHDU().write(1, dim, array);
}

void write_pimg(string filename, int x_dim, int y_dim, gsl_vector *data)
{

  //see above


  long naxis = 2;
  long naxes[2] = {x_dim, y_dim};
  int dim = x_dim * y_dim;
  filename = "!"+filename;

  std::auto_ptr<FITS> img(0);
  img.reset(new FITS(filename,DOUBLE_IMG,naxis,naxes));


  std::valarray<double> array(dim); 
  for (int i = 0; i < dim; i++)
    {
      array[i] = gsl_vector_get(data,i);
    }
  
  
  img->pHDU().write(1, dim, array);
}


void write_pimgint(string filename, int x_dim, int y_dim, gsl_vector_int *data)
{

  long naxis = 2;
  long naxes[2] = {x_dim, y_dim};
  int dim = x_dim * y_dim;
  filename = "!"+filename;

  std::auto_ptr<FITS> img(0);
  img.reset(new FITS(filename,32,naxis,naxes));

  std::valarray<int> array(dim); 
  for (int i = 0; i < dim; i++)
    {
      array[i] = gsl_vector_int_get(data,i);
    }
  
    
  img->pHDU().write(1, dim, array);
}

void write_imge(string filename, string extname, gsl_matrix *data)
{

  int x_dim = data->size2;
  int y_dim = data->size1; 
  int dim = x_dim * y_dim;

  std::auto_ptr<FITS> img(0);
  img.reset(new FITS(filename,Write));


  std::valarray<double> array(dim); 
  for (int i = 0; i < y_dim; i++)
    {
      for(int j = 0; j < x_dim; j++)
	{
	  array[i*x_dim + j] = gsl_matrix_get(data, i, j);
	}
    }
  
  std::vector<long> extAx(2, 0);
  extAx.front() = x_dim;
  extAx.back() = y_dim;
  ExtHDU* imageExt = img->addImage(extname, DOUBLE_IMG, extAx);
  imageExt->write(1, dim, array);
}

void write_imgeint(string filename, string extname, gsl_matrix_int *data)
{
  int x_dim = data->size2;
  int y_dim = data->size1;
  int dim = x_dim * y_dim;
  std::auto_ptr<FITS> img(0);
  img.reset(new FITS(filename,Write));
  
  std::valarray<int> array(dim); 
  for (int i = 0; i < y_dim; i++)
    {
      for(int j = 0; j < x_dim; j++)
	{
	  array[i*x_dim + j] = gsl_matrix_int_get(data, i, j);
	}
    }
  
  std::vector<long> extAx(2, 0);
  extAx.front() = x_dim;
  extAx.back() = y_dim;
  ExtHDU* imageExt = img->addImage(extname, 32, extAx);
  imageExt->write(1, dim, array);
}

void write_imge(string filename, string extname, int x_dim, int y_dim, gsl_vector *data)
{
  int dim = x_dim * y_dim;
  std::auto_ptr<FITS> img(0);
  img.reset(new FITS(filename,Write));
  
  std::valarray<double> array(dim); 
  for (int i = 0; i < dim; i++)
    {
      array[i] = gsl_vector_get(data,i);
    }
  
  std::vector<long> extAx(2, 0);
  extAx.front() = x_dim;
  extAx.back() = y_dim;
  ExtHDU* imageExt = img->addImage(extname, DOUBLE_IMG, extAx);
  imageExt->write(1, dim, array);
}

void write_imgeint(string filename, string extname, int x_dim, int y_dim, gsl_vector_int *data)
{
  int dim = x_dim * y_dim;
  std::auto_ptr<FITS> img(0);
  img.reset(new FITS(filename,Write));
  
  std::valarray<int> array(dim); 
  for (int i = 0; i < dim; i++)
    {
      array[i] = gsl_vector_int_get(data,i);
    }
  
  std::vector<long> extAx(2, 0);
  extAx.front() = x_dim;
  extAx.back() = y_dim;
  ExtHDU* imageExt = img->addImage(extname, 32, extAx);
  imageExt->write(1, dim, array);
}

void read_pimg(string filename, gsl_matrix *data)
{
  //reads the double pHDU of a fits file

  std::auto_ptr<FITS> pInfile(new FITS(filename, Read, true));
  PHDU& img = pInfile->pHDU();

  
  std::valarray<double> contents;
  
  
  img.readAllKeys();
  img.read(contents);
  int ax1(img.axis(0));
  int ax2(img.axis(1));



  for (int i = 0; i < ax2; i++)
    {
      for(int j = 0; j < ax1; j++)
	{
	  gsl_matrix_set(data, i, j, contents[i*ax1 + j]);
	}
    }
  
}

void read_pimgint(string filename, gsl_matrix_int *data)
{

  //same in integer

  std::auto_ptr<FITS> pInfile(new FITS(filename, Read, true));
  PHDU& img = pInfile->pHDU();

 
  std::valarray<int> contents;


  img.readAllKeys();
  img.read(contents);
  int ax1(img.axis(0));
  int ax2(img.axis(1));

  
  
  for (int i = 0; i < ax2; i++)
    {
      for(int j = 0; j < ax1; j++)
	{
	  gsl_matrix_int_set(data, i, j, contents[i*ax1 + j]);
	}
    }
  
}

void read_pimg(string filename, gsl_vector *data)
{
  //see above

  std::auto_ptr<FITS> pInfile(new FITS(filename, Read, true));
  PHDU& img = pInfile->pHDU();
  
  
  std::valarray<double> contents;
  

  img.readAllKeys();
  img.read(contents);
  int ax1(img.axis(0));
  int ax2(img.axis(1));
  int dim = ax1 * ax2;



  for (int i = 0; i < dim; i++)
    {
      gsl_vector_set(data, i, contents[i]);
    }

}

void read_pimgint(string filename, gsl_vector_int *data)
{

  //see above


  std::auto_ptr<FITS> pInfile(new FITS(filename, Read, true));
  PHDU& img = pInfile->pHDU();

 
  std::valarray<int> contents;


  img.readAllKeys();
  img.read(contents);
  int ax1(img.axis(0));
  int ax2(img.axis(1));
  int dim = ax1 * ax2;
  


  for (int i = 0; i < dim; i++)
    {
      gsl_vector_int_set(data, i, contents[i]);
    }

}

void read_imge(string filename, string extname, gsl_matrix *data)
{

  //reads a double image extension of a fits file

  std::auto_ptr<FITS> pInfile(new FITS(filename, Read, true));
  ExtHDU& img = pInfile->extension(extname);
  std::valarray<double> contents;
  img.readAllKeys();
  img.read(contents);
  int ax1(img.axis(0));
  int ax2(img.axis(1));

  for (int i = 0; i < ax2; i++)
    {
      for(int j = 0; j < ax1; j++)
	{
	  gsl_matrix_set(data, i, j, contents[i*ax1 + j]);
	}
    }
}

void read_imgeint(string filename, string extname, gsl_matrix_int *data)
{

  //same in integer

  std::auto_ptr<FITS> pInfile(new FITS(filename, Read, true));
  ExtHDU& img = pInfile->extension(extname);
  std::valarray<int> contents;
  img.readAllKeys();
  img.read(contents);
  int ax1(img.axis(0));
  int ax2(img.axis(1));

  for (int i = 0; i < ax2; i++)
    {
      for(int j = 0; j < ax1; j++)
	{
	  gsl_matrix_int_set(data, i, j, contents[i*ax1 + j]);
	}
    }
}

void read_imge(string filename, string extname, gsl_vector *data)
{

  //see above
  std::auto_ptr<FITS> pInfile(new FITS(filename, Read, true));
  ExtHDU& img = pInfile->extension(extname);
  std::valarray<double> contents;
  img.readAllKeys();
  img.read(contents);
  int ax1(img.axis(0));
  int ax2(img.axis(1));
  int dim = ax1 * ax2;

  for (int i = 0; i < dim; i++)
    {
      gsl_vector_set(data,i, contents[i]);
    }
}

void read_imgeint(string filename, string extname, gsl_vector_int *data)
{

  //see above

  std::auto_ptr<FITS> pInfile(new FITS(filename, Read, true));
  ExtHDU& img = pInfile->extension(extname);
  std::valarray<int> contents;
  img.readAllKeys();
  img.read(contents);
  int ax1(img.axis(0));
  int ax2(img.axis(1));
  int dim = ax1 * ax2;
  for (int i = 0; i < dim; i++)
    {
      gsl_vector_int_set(data,i, contents[i]);
    }
}

void write_header(string filename,string Key, double value, string description)
{

  std::auto_ptr<FITS> img(0);
  img.reset(new FITS(filename,Write));

  img->pHDU().addKey(Key,value,description);

}

void write_header(string filename,string Key, int value, string description)
{

  std::auto_ptr<FITS> img(0);
  img.reset(new FITS(filename,Write));

  img->pHDU().addKey(Key,value,description);

}

void write_header(string filename,string Key, string value, string description)
{

  std::auto_ptr<FITS> img(0);
  img.reset(new FITS(filename,Write));

  img->pHDU().addKey(Key,value,description);

}

double read_doubleheader(string filename, string name)
{
  double read;

  std::auto_ptr<FITS> pInfile(new FITS(filename, Read, true));
  PHDU& img = pInfile->pHDU();

  img.readKey(name,read);

  return read;
}

int read_intheader(string filename, string name)
{
  int read;

  std::auto_ptr<FITS> pInfile(new FITS(filename, Read, true));
  PHDU& img = pInfile->pHDU();

  img.readKey(name,read);

  return read;
}

string read_stringheader(string filename, string name)
{
  string read;

  std::auto_ptr<FITS> pInfile(new FITS(filename, Read, true));
  PHDU& img = pInfile->pHDU();

  img.readKey(name,read);

  return read;
}











