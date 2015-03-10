/***
    NON-public source file.
    Author: Julian Merten
    Inst.: ITA/ZAH Heidelberg
    2010
***/

#include <cat_reader.h>


CatReader::CatReader(const string &filename)
{

    int row_counter, col_counter;
    vector<double> col_check;
    string line;

    ifstream input(filename.c_str());

    if(!input.is_open())
	{
	    throw invalid_argument("Initial CatReader input invalid");
	}
    row_counter = 0;
    col_counter = 0;
    while(getline(input,line))
	{
	    col_check.clear();
	    if(line.find_first_of("01234567890.+-e",line.find_first_not_of(" ")) < line.find_first_not_of("01234567890.+-e",line.find_first_not_of(" ")))
		{
		    read_doubles(line,col_check);
		    if(col_check.size() > col_counter)
			{
			    col_counter = col_check.size();
			}
		    if(col_check.size() > 0)
			{
			    row_counter++;
			}
		}
	}
    input.close();

    data = gsl_matrix_calloc(row_counter,col_counter);

    input.open(filename.c_str());

    if(!input)
	{
	    throw invalid_argument("Initial CatReader input invalid");
	}

    row_counter = 0;

    while(getline(input,line))
	{
	    col_check.clear();
	    if(line.find_first_of("01234567890.+-e",line.find_first_not_of(" ")) < line.find_first_not_of("01234567890.+-e",line.find_first_not_of(" ")))
		{
		    read_doubles(line,col_check);
		    if(col_check.size() != data->size2)
			{
			    cout <<"Warning: There seems to be an invalid line in your CatReader input." <<endl;
			}
		    for(int i = 0; i < col_check.size(); i++)
			{
			    gsl_matrix_set(data,row_counter,i,col_check[i]);
			}
		    row_counter++;
		}
	}
    cat_counter = 1;
    input.close();
    col_check.clear();
    col_names.resize(data->size2, " ");
}

CatReader::CatReader(const string &filename, string symbol)
{
    int row_counter, col_counter;
    vector<double> col_check;
    string line, subline;
    col_names.clear();

    ifstream input(filename.c_str());

    if(!input.is_open())
	{
	    throw invalid_argument("Initial CatReader input invalid");
	}
    row_counter = 0;
    col_counter = 0;
    while(getline(input,line))
	{
	    col_check.clear();
	    if(line.find_first_of("01234567890.+-e",line.find_first_not_of(" ")) < line.find_first_not_of("01234567890.+-e",line.find_first_not_of(" ")))
		{
		    read_doubles(line,col_check);
		    if(col_check.size() > col_counter)
			{
			    col_counter = col_check.size();
			}
		    if(col_check.size() > 0)
			{
			    row_counter++;
			}
		}
	}
    input.close();

    data = gsl_matrix_calloc(row_counter,col_counter);

    input.open(filename.c_str());

    if(!input)
	{
	    throw invalid_argument("Initial CatReader input invalid");
	}

    row_counter = 0;

    while(getline(input,line))
	{
	    col_check.clear();
	    if(line.find_first_of("01234567890.+-e",line.find_first_not_of(" ")) < line.find_first_not_of("01234567890.+-e",line.find_first_not_of(" ")))
		{
		    read_doubles(line,col_check);
		    if(col_check.size() != data->size2)
			{
			    cout <<"Warning: There seems to be an invalid line in your CatReader input." <<endl;
			}
		    for(int i = 0; i < col_check.size(); i++)
			{
			    gsl_matrix_set(data,row_counter,i,col_check[i]);
			}
		    row_counter++;
		}
	    else
		{
		    while(line.find(symbol,0) != string::npos && col_names.size() <= data->size2)
			{
			    line.erase(0,line.find(symbol,0));
			    subline = line.substr(line.find_first_not_of(symbol+" "),line.find_first_of(" "+symbol,line.find_first_not_of(symbol+" ")));
			    line.erase(0,subline.size());
			    col_names.push_back(subline);
			}
		}
	}
    cat_counter = 1;
    input.close();
    col_check.clear();
    col_names.resize(data->size2, " ");
}


CatReader::CatReader(const string &filename, int row_limit, int col_limit)
{

    int row_counter;
    vector<double> col_check;
    string line;

    ifstream input(filename.c_str());

    if(!input)
	{
	    throw invalid_argument("Initial CatReader input invalid");
	}
    row_counter = 0;
    while(getline(input,line))
	{
	    col_check.clear();
	    if(line.find_first_of("01234567890.+-e",line.find_first_not_of(" ")) < line.find_first_not_of("01234567890.+-e",line.find_first_not_of(" ")) && row_counter < row_limit)
		{
		    read_doubles(line, col_limit, col_check);
		    if(col_check.size() > 0)
			{
			    row_counter++;
			}
		}
	}
    input.close();

    data = gsl_matrix_calloc(row_counter,col_limit);

    input.open(filename.c_str());

    if(!input)
	{
	    throw invalid_argument("Initial CatReader input invalid");
	}

    row_counter = 0;

    while(getline(input,line))
	{
	    col_check.clear();
	    if(line.find_first_of("01234567890.+-e",line.find_first_not_of(" ")) < line.find_first_not_of("01234567890.+-e",line.find_first_not_of(" ")) && row_counter < row_limit)
		{
		    read_doubles(line,col_limit, col_check);
		    if(col_check.size() != data->size2)
			{
			    cout <<"Warning: There seems to be an invalid line in your CatReader input." <<endl;
			}
		    for(int i = 0; i < col_check.size(); i++)
			{
			    gsl_matrix_set(data,row_counter,i,col_check[i]);
			}
		    row_counter++;
		}
	}
    cat_counter = 1;
    input.close();
    col_check.clear();
    col_names.resize(data->size2, " ");
}

CatReader::CatReader(const string &filename, int row_limit, int col_limit, string symbol)
{

    int row_counter;
    vector<double> col_check;
    string line, subline;
    col_names.clear();

    ifstream input(filename.c_str());

    if(!input)
	{
	    throw invalid_argument("Initial CatReader input invalid");
	}
    row_counter = 0;
    while(getline(input,line))
	{
	    col_check.clear();
	    if(line.find_first_of("01234567890.+-e",line.find_first_not_of(" ")) < line.find_first_not_of("01234567890.+-e",line.find_first_not_of(" ")) && row_counter < row_limit)
		{
		    read_doubles(line, col_limit, col_check);
		    if(col_check.size() > 0)
			{
			    row_counter++;
			}
		}
	}
    input.close();

    data = gsl_matrix_calloc(row_counter,col_limit);

    input.open(filename.c_str());

    if(!input)
	{
	    throw invalid_argument("Initial CatReader input invalid");
	}

    row_counter = 0;

    while(getline(input,line))
	{
	    col_check.clear();
	    if(line.find_first_of("01234567890.+-e",line.find_first_not_of(" ")) < line.find_first_not_of("01234567890.+-e",line.find_first_not_of(" ")) && row_counter < row_limit)
		{
		    read_doubles(line,col_limit, col_check);
		    if(col_check.size() != data->size2)
			{
			    cout <<"Warning: There seems to be an invalid line in your CatReader input." <<endl;
			}
		    for(int i = 0; i < col_check.size(); i++)
			{
			    gsl_matrix_set(data,row_counter,i,col_check[i]);
			}
		    row_counter++;
		}
	    else
		{
		    while(line.find(symbol,0) != string::npos && col_names.size() <= data->size2)
			{
			    line.erase(0,line.find(symbol,0));
			    subline = line.substr(line.find_first_not_of(symbol+" "),line.find_first_of(" "+symbol,line.find_first_not_of(symbol+" ")));
			    line.erase(0,subline.size());
			    col_names.push_back(subline);
			}
		}
	}
    cat_counter = 1;
    input.close();
    col_check.clear();
    col_names.resize(data->size2, " ");
}

CatReader::~CatReader()
{

    col_names.clear();
    gsl_matrix_free(data);

}

int CatReader::show_rowdim()
{
    return data->size1;
}

int CatReader::show_coldim()
{

    return data->size2;
}


void CatReader::read(const string &filename)
{

    int row_counter, col_counter;
    vector<double> col_check;
    string line;

    ifstream input(filename.c_str());

    if(!input.is_open())
	{
	    throw invalid_argument("Additional CatReader input invalid");
	}
    row_counter = 0;
    col_counter = 0;
    while(getline(input,line))
	{
	    col_check.clear();
	    if(line.find_first_of("01234567890.+-e",line.find_first_not_of(" ")) < line.find_first_not_of("01234567890.+-e",line.find_first_not_of(" ")))
		{
		    read_doubles(line,col_check);
		    if(col_check.size() > col_counter)
			{
			    col_counter = col_check.size();
			}
		    if(col_check.size() > 0)
			{
			    row_counter++;
			}
		}
	}
    input.close();

    if(row_counter != data->size1)
	{
	    throw invalid_argument("Row number of new cat is invalid");
	}

    gsl_matrix *newdata = gsl_matrix_calloc(row_counter,col_counter);

    input.open(filename.c_str());

    if(!input.is_open())
	{
	    throw invalid_argument("Additional CatReader input invalid");
	}

    row_counter = 0;

    while(getline(input,line))
	{
	    col_check.clear();
	    if(line.find_first_of("01234567890.+-e",line.find_first_not_of(" ")) < line.find_first_not_of("01234567890.+-e",line.find_first_not_of(" ")))
		{
		    read_doubles(line,col_check);
		    if(col_check.size() != newdata->size2)
			{
			    cout <<"Warning: There seems to be an invalid line in your CatReader input." <<endl;
			}
		    for(int i = 0; i < col_check.size(); i++)
			{
			    gsl_matrix_set(newdata,row_counter,i,col_check[i]);
			}
		    row_counter++;
		}
	}

    gsl_matrix *olddata = gsl_matrix_calloc(data->size1,data->size2);
    gsl_matrix_memcpy(olddata,data);
    gsl_matrix_free(data);
    data = gsl_matrix_calloc(row_counter,olddata->size2+newdata->size2);

    for(int i = 0; i < data->size1; i++)
	{
	    for(int j = 0; j < data->size2; j++)
		{
		    if(j < olddata->size2)
			{
			    gsl_matrix_set(data,i,j,gsl_matrix_get(olddata,i,j));
			}
		    else
			{
			    gsl_matrix_set(data,i,j,gsl_matrix_get(newdata,i,j-olddata->size2));
			}
		}
	}

    gsl_matrix_free(newdata);
    gsl_matrix_free(olddata);
    cat_counter++;
    col_names.resize(data->size2, " ");

}

void CatReader::read(const string &filename, string symbol)
{

    int row_counter, col_counter;
    vector<double> col_check;
    string line, subline;

    ifstream input(filename.c_str());

    if(!input.is_open())
	{
	    throw invalid_argument("Additional CatReader input invalid");
	}
    row_counter = 0;
    col_counter = 0;
    while(getline(input,line))
	{
	    col_check.clear();
	    if(line.find_first_of("01234567890.+-e",line.find_first_not_of(" ")) < line.find_first_not_of("01234567890.+-e",line.find_first_not_of(" ")))
		{
		    read_doubles(line,col_check);
		    if(col_check.size() > col_counter)
			{
			    col_counter = col_check.size();
			}
		    if(col_check.size() > 0)
			{
			    row_counter++;
			}
		}
	}
    input.close();

    if(row_counter != data->size1)
	{
	    throw invalid_argument("Row number of new cat is invalid");
	}

    gsl_matrix *newdata = gsl_matrix_calloc(row_counter,col_counter);

    input.open(filename.c_str());

    if(!input.is_open())
	{
	    throw invalid_argument("Additional CatReader input invalid");
	}

    row_counter = 0;

    while(getline(input,line))
	{
	    col_check.clear();
	    if(line.find_first_of("01234567890.+-e",line.find_first_not_of(" ")) < line.find_first_not_of("01234567890.+-e",line.find_first_not_of(" ")))
		{
		    read_doubles(line,col_check);
		    if(col_check.size() != newdata->size2)
			{
			    cout <<"Warning: There seems to be an invalid line in your CatReader input." <<endl;
			}
		    for(int i = 0; i < col_check.size(); i++)
			{
			    gsl_matrix_set(newdata,row_counter,i,col_check[i]);
			}
		    row_counter++;
		}
	    else
		{
		    while(line.find(symbol,0) != string::npos && col_names.size() <= data->size2+newdata->size2)
			{
			    line.erase(0,line.find(symbol,0));
			    subline = line.substr(line.find_first_not_of(symbol+" "),line.find_first_of(" "+symbol,line.find_first_not_of(symbol+" ")));
			    line.erase(0,subline.size());
			    col_names.push_back(subline);
			}
		}
	}

    gsl_matrix *olddata = gsl_matrix_calloc(data->size1,data->size2);
    gsl_matrix_memcpy(olddata,data);
    gsl_matrix_free(data);
    data = gsl_matrix_calloc(row_counter,olddata->size2+newdata->size2);

    for(int i = 0; i < data->size1; i++)
	{
	    for(int j = 0; j < data->size2; j++)
		{
		    if(j < olddata->size2)
			{
			    gsl_matrix_set(data,i,j,gsl_matrix_get(olddata,i,j));
			}
		    else
			{
			    gsl_matrix_set(data,i,j,gsl_matrix_get(newdata,i,j-olddata->size2));
			}
		}
	}

    gsl_matrix_free(newdata);
    gsl_matrix_free(olddata);
    cat_counter++;
    col_names.resize(data->size2, " ");

}


void CatReader::read(const string &filename, int row_limit, int col_limit)
{

    int row_counter;
    vector<double> col_check;
    string line;

    ifstream input(filename.c_str());

    if(!input.is_open())
	{
	    throw invalid_argument("Additional CatReader input invalid");
	}
    row_counter = 0;
    while(getline(input,line))
	{
	    col_check.clear();
	    if(line.find_first_of("01234567890.+-e",line.find_first_not_of(" ")) < line.find_first_not_of("01234567890.+-e",line.find_first_not_of(" ")) && row_counter < row_limit)
		{
		    read_doubles(line, col_limit, col_check);
		    if(col_check.size() > 0)
			{
			    row_counter++;
			}
		}
	}
    input.close();

    if(row_counter != data->size1)
	{
	    throw invalid_argument("Row number of new cat is invalid");
	}

    gsl_matrix *newdata = gsl_matrix_calloc(row_counter,col_limit);

    input.open(filename.c_str());

    if(!input.is_open())
	{
	    throw invalid_argument("Additional CatReader input invalid");
	}

    row_counter = 0;

    while(getline(input,line))
	{
	    col_check.clear();
	    if(line.find_first_of("01234567890.+-e",line.find_first_not_of(" ")) < line.find_first_not_of("01234567890.+-e",line.find_first_not_of(" ")) && row_counter < row_limit)
		{
		    read_doubles(line,col_limit, col_check);
		    if(col_check.size() != newdata->size2)
			{
			    cout <<"Warning: There seems to be an invalid line in your CatReader input." <<endl;
			}
		    for(int i = 0; i < col_check.size(); i++)
			{
			    gsl_matrix_set(newdata,row_counter,i,col_check[i]);
			}
		    row_counter++;
		}
	}

    input.close();

    gsl_matrix *olddata = gsl_matrix_calloc(data->size1,data->size2);
    gsl_matrix_memcpy(olddata,data);
    gsl_matrix_free(data);
    data = gsl_matrix_calloc(row_counter,olddata->size2+newdata->size2);

    for(int i = 0; i < data->size1; i++)
	{
	    for(int j = 0; j < data->size2; j++)
		{
		    if(j < olddata->size2)
			{
			    gsl_matrix_set(data,i,j,gsl_matrix_get(olddata,i,j));
			}
		    else
			{
			    gsl_matrix_set(data,i,j,gsl_matrix_get(newdata,i,j-olddata->size2));
			}
		}
	}

    gsl_matrix_free(newdata);
    gsl_matrix_free(olddata);
    cat_counter++;
    col_names.resize(data->size2," ");
}

void CatReader::read(const string &filename, int row_limit, int col_limit, string symbol)
{

    int row_counter;
    vector<double> col_check;
    string line, subline;

    ifstream input(filename.c_str());

    if(!input.is_open())
	{
	    throw invalid_argument("Additional CatReader input invalid");
	}
    row_counter = 0;
    while(getline(input,line))
	{
	    col_check.clear();
	    if(line.find_first_of("01234567890.+-e",line.find_first_not_of(" ")) < line.find_first_not_of("01234567890.+-e",line.find_first_not_of(" ")) && row_counter < row_limit)
		{
		    read_doubles(line, col_limit, col_check);
		    if(col_check.size() > 0)
			{
			    row_counter++;
			}
		}
	}
    input.close();

    if(row_counter != data->size1)
	{
	    throw invalid_argument("Row number of new cat is invalid");
	}

    gsl_matrix *newdata = gsl_matrix_calloc(row_counter,col_limit);

    input.open(filename.c_str());

    if(!input.is_open())
	{
	    throw invalid_argument("Additional CatReader input invalid");
	}

    row_counter = 0;

    while(getline(input,line))
	{
	    col_check.clear();
	    if(line.find_first_of("01234567890.+-e",line.find_first_not_of(" ")) < line.find_first_not_of("01234567890.+-e",line.find_first_not_of(" ")) && row_counter < row_limit)
		{
		    read_doubles(line,col_limit, col_check);
		    if(col_check.size() != newdata->size2)
			{
			    cout <<"Warning: There seems to be an invalid line in your CatReader input." <<endl;
			}
		    for(int i = 0; i < col_check.size(); i++)
			{
			    gsl_matrix_set(newdata,row_counter,i,col_check[i]);
			}
		    row_counter++;
		}
	    else
		{
		    while(line.find(symbol,0) != string::npos && col_names.size() <= data->size2+newdata->size2)
			{
			    line.erase(0,line.find(symbol,0));
			    subline = line.substr(line.find_first_not_of(symbol+" "),line.find_first_of(" "+symbol,line.find_first_not_of(symbol+" ")));
			    line.erase(0,subline.size());
			    col_names.push_back(subline);
			}
		}
	}

    input.close();

    gsl_matrix *olddata = gsl_matrix_calloc(data->size1,data->size2);
    gsl_matrix_memcpy(olddata,data);
    gsl_matrix_free(data);
    data = gsl_matrix_calloc(row_counter,olddata->size2+newdata->size2);

    for(int i = 0; i < data->size1; i++)
	{
	    for(int j = 0; j < data->size2; j++)
		{
		    if(j < olddata->size2)
			{
			    gsl_matrix_set(data,i,j,gsl_matrix_get(olddata,i,j));
			}
		    else
			{
			    gsl_matrix_set(data,i,j,gsl_matrix_get(newdata,i,j-olddata->size2));
			}
		}
	}

    gsl_matrix_free(newdata);
    gsl_matrix_free(olddata);
    cat_counter++;
    col_names.resize(data->size2," ");
}


void CatReader::set_names(int col, string name)
{

    if(col < col_names.size())
	{
	    col_names[col] = name;
	}
    else
	{
	    throw invalid_argument("Invalid column to a name to");
	}
}

string CatReader::show_names(int col)
{

    if(col < col_names.size())
	{
	    return col_names[col];
	}
    else
	{
	    throw invalid_argument("Invalid column to a name to");
	}
}

int CatReader::show_counter()
{

    return cat_counter;

}

gsl_matrix* CatReader::show_data()
{

    return data;
}

void CatReader::write_data(gsl_matrix *output)
{

    if(output->size1 == data->size1 && output->size2 == data->size2)
	{
	    gsl_matrix_memcpy(output,data);
	}

    else
	{
	    throw invalid_argument("Matrix dims must match to write cat");
	}

}














