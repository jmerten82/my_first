/***
    NON-public source file.
    Author: Julian Merten
    Inst.: ITA/ZAH Heidelberg
    2010
***/

#include <cat_writer.h>

CatWriter::CatWriter(CatReader *cat_input)
{

    input = cat_input;
    row_flags.resize(input->show_rowdim(),true);

    for(int i = 0; i < input->show_coldim(); i++)
	{
	    col_flags.push_back(i);
	    names.push_back(input->show_names(i));
	}
    num_rows = row_flags.size();
    num_columns = col_flags.size();
    sep = "\t";
    mode = "scientific";
    write_info = false;
    alloc = false;
    update_needed = true;
    WCS = false;
}

CatWriter::CatWriter(CatReader *cat_input, string col_select)
{

    input = cat_input;
    row_flags.resize(input->show_rowdim(),true);

    read_sequence(col_select,col_flags);
    col_flags.resize(input->show_coldim(),-1);

    num_rows = row_flags.size();
    num_columns = 0;

    for(int i = 0; i < col_flags.size(); i++)
	{
	    if(col_flags[i] != -1)
		{
		    num_columns++;
		    names.push_back(input->show_names(col_flags[i]));
		}
	}

    sep = "\t";
    mode = "scientific";
    write_info = false;
    alloc = false;
    update_needed = true;
    WCS = false;
}

CatWriter::CatWriter(CatReader *cat_input, string col_select, string row_select)
{

    input = cat_input;

    row_flags.resize(input->show_rowdim(),false);

    vector<int> helper;
    read_sequence(row_select,helper);

    for(int i = 0; i < helper.size(); i++)
	{
	    row_flags[helper[i]] = true;
	}
    helper.clear();

    if(col_select == "all")
	{
	    col_flags.clear();
	    for(int i = 0; i < input->show_coldim(); i++)
		{
		    col_flags.push_back(i);
		}
	}
    else
	{
	    col_flags.resize(input->show_coldim(),-1);
	    read_sequence(col_select,col_flags);
	}

    num_rows = 0;
    num_columns = 0;

    for(int i = 0; i < row_flags.size(); i++)
	{
	    if(row_flags[i])
		{
		    num_rows++;
		}
	}
    for(int i = 0; i < col_flags.size(); i++)
	{
	    if(col_flags[i] != -1)
		{
		    num_columns++;
		    names.push_back(input->show_names(col_flags[i]));
		}
	}

    sep = "\t";
    mode = "scientific";
    write_info = false;
    alloc = false;
    update_needed = true;
    WCS = false;
}

CatWriter::~CatWriter()
{

    row_flags.clear();
    col_flags.clear();
    names.clear();

    if(alloc)
	{
	    gsl_matrix_free(data);
	}
}

void CatWriter::set_col(int col, int value)
{

    if(col >= col_flags.size() || value >= col_flags.size())
	{
	    throw invalid_argument("Invalid column selection");
	}

    col_flags[col] = value;
    num_columns = 0;
    names.clear();
    for(int i = 0; i < col_flags.size(); i++)
	{
	    if(col_flags[i] != -1)
		{
		    num_columns++;
		    names.push_back(input->show_names(col_flags[i]));
		}
	}
    update_needed = true;
}
void CatWriter::set_col(string sequence)
{
    col_flags.clear();
    read_sequence(sequence,col_flags);
    col_flags.resize(input->show_coldim(),-1);
    names.clear();

    num_columns = 0;
    for(int i = 0; i < col_flags.size(); i++)
	{
	    if(col_flags[i] != -1)
		{
		    num_columns++;
		    names.push_back(input->show_names(col_flags[i]));
		}
	}
    update_needed = true;
}

int CatWriter::show_col(int col)
{

    if(col >= col_flags.size())
	{
	    throw invalid_argument("Invalid column selection");
	}
    return col_flags[col];
}

void CatWriter::set_row(int row, bool value)
{

    if(row >= row_flags.size() || value >= row_flags.size())
	{
	    throw invalid_argument("Invalid row selection");
	}

    row_flags[row] = value;
    num_rows = 0;
    for(int i = 0; i < row_flags.size(); i++)
	{
	    if(row_flags[i])
		{
		    num_rows++;
		}
	}
    update_needed = true;
}
void CatWriter::set_row(string sequence)
{
    row_flags.clear();
    row_flags.resize(input->show_coldim(),false);
    vector<int> helper;
    read_sequence(sequence,helper);
    for(int i = 0; i < helper.size(); i++)
	{
	    row_flags[helper[i]] = true;
	}

    num_rows = 0;
    for(int i = 0; i < row_flags.size(); i++)
	{
	    if(row_flags[i])
		{
		    num_rows++;
		}
	}
    update_needed = true;
}

bool CatWriter::show_row(int row)
{

    if(row >= row_flags.size())
	{
	    throw invalid_argument("Invalid row selection");
	}
    return row_flags[row];
}

void CatWriter::set_mode(string symbol, string selection)
{

    if(symbol.size() > 0)
	{
	    sep = symbol;
	}

    if(selection == "normal" || selection == "scientific" || selection == "int")
	{
	    mode = selection;
	}
    else
	{
	    throw invalid_argument("Invalid selection for set_mode");
	}
}

int CatWriter::show_rowdim()
{
    return num_rows;
}

int CatWriter::show_coldim()
{
    return num_columns;
}

void CatWriter::set_infomode(bool value)
{

    write_info = value;

}

string CatWriter::show_mode()
{
    return mode;
}

bool CatWriter::show_infomode()
{

    return write_info;
}

void CatWriter::write(const string &filename)
{

    ofstream output(filename.c_str());

    if(write_info)
	{
	    for(int i = 0; i < num_columns; i++)
		{
		    if(i != num_columns-1)
			{
			    output <<"#"<<i <<" " <<names[i] <<"\t" <<flush;
			}
		    else
			{

			    output <<"#"<<i <<" " <<names[i] <<endl;
			}
		}
	}



    if(mode == "scientific")
	{
	    output <<scientific;
	    output <<setprecision(7);
	    output <<showpoint;
	}
    else if(mode == "normal")
	{
	    output <<fixed;
	}


    if(!alloc)
	{
	    
	    for(int i = 0; i < num_rows; i++)
		{
		    if(row_flags[i])
			{
			    for(int j = 0; j < num_columns; j++)
				{
				    if(col_flags[j] != -1)
					{
					    if(j != 0)
						{
						    output <<sep <<gsl_matrix_get(input->show_data(),i,col_flags[j]) <<flush;
						}
					    else
						{
						    output <<gsl_matrix_get(input->show_data(),i,col_flags[j]) <<flush;
						}
					}
				}
			    output <<endl;
			}
		}
	}
    else if(update_needed)
	{
	    create_copy();
	    for(int i = 0; i < data->size1; i++)
		{
		    for(int j = 0; j < data->size2; j++)
			{
			    if(j != 0)
				{
				    output <<sep <<gsl_matrix_get(data,i,j) <<flush;
				}
			    else
				{
				    output <<gsl_matrix_get(data,i,j) <<flush;
				}
			}
		    output <<endl;
		}
	}
    else
	{
	    for(int i = 0; i < data->size1; i++)
		{
		    for(int j = 0; j < data->size2; j++)
			{
			    if(j != 0)
				{
				    output <<sep <<gsl_matrix_get(data,i,j) <<flush;
				}
			    else
				{
				    output <<gsl_matrix_get(data,i,j) <<flush;
				}
			}
		    output <<endl;
		}
	}

			    

}

void CatWriter::write(gsl_matrix *output)
{

    if(output->size1 != num_rows && output->size2 != num_columns)
	{
	    throw invalid_argument("Matrix dims must match for cat_writer write");
	}

    if(!alloc)
	{
	    
	    for(int i = 0; i < row_flags.size(); i++)
		{
		    if(row_flags[i])
			{
			    for(int j = 0; j < col_flags.size(); j++)
				{
				    if(col_flags[j] != -1)
					{
					    gsl_matrix_set(output,i,j,gsl_matrix_get(input->show_data(),i,col_flags[j]));
					}
				}
			}
		}
	}
    else if(update_needed)
	{
	    create_copy();
	    gsl_matrix_memcpy(output,data);
	}
    else
	{
	    gsl_matrix_memcpy(output,data);
	}
}

void CatWriter::create_copy()
{

    if(alloc)
	{
	    gsl_matrix_free(data);
	}

    data = gsl_matrix_calloc(num_rows,num_columns);


    for(int i = 0; i < row_flags.size(); i++)
	{
	    if(row_flags[i])
		{
		    for(int j = 0; j < col_flags.size(); j++)
			{
			    if(col_flags[j] != -1)
				{
				    gsl_matrix_set(data,i,j,gsl_matrix_get(input->show_data(),i,col_flags[j]));
				}
			}
		}
	}
    alloc = true;
    update_needed = false;
}

gsl_matrix* CatWriter::show_data()
{

    if(!alloc || update_needed)
	{
	    create_copy();
	}
    else
	{
	    return data;
	}
}

void CatWriter::set_WCS(bool flag)
{
    WCS = flag;
}

bool CatWriter::show_WCS()
{
    return WCS;
}

void CatWriter::write_DS9_points(const string &filename, int col1, int col2, string selection)
{

    if(!alloc || update_needed)
	{
	    create_copy();

     	}

    if(col1 >= data->size2 || col2 >= data->size2)
	{
	    throw invalid_argument("Invalid column selection for DS9 write");
	}

    ofstream output(filename.c_str());

    output <<"# Region file format: DS9 version 4.0" <<endl;
    output <<"# SaWLens catalogue writer: ASCII->DS9 conversion." <<endl;

    string mode;

    if(selection == "circle")
	{
	    mode = "circle";
	}
    else if(selection == "box")
	{
	    mode = "box";
	}
    else if(selection == "diamond")
	{
	    mode = "diamond";
	}
    else if(selection == "cross")
	{
	    mode = "cross";
	}
    else if(selection == "xpoint")
	{
	    mode = "x";
	}
    else if(selection == "boxcircle")
	{
	    mode = "boxcircle";
	}
    else
	{
	    throw invalid_argument("Invalid point selection in DS9 write");
	}

    output <<"global color=white dashlist=8 3 width=1 font=\"helvetica 10 normal\" select=1 highlite=0 dash=0 fixed=1 edit=0 move=0 delete=1 include=1 source=0" <<endl;

    if(WCS)
	{
	    output <<"fk5" <<endl;
	}
    else
	{
	    output <<"linear" <<endl;
	}


    for(int i = 0; i < data->size1; i++)
	{
	    output <<"point(" <<gsl_matrix_get(data,i,col1) <<"," <<gsl_matrix_get(data,i,col2) <<") # point=" <<mode <<" tag={WL}" <<endl;
	}

    output.close();
}

void CatWriter::write_DS9_circles(const string &filename, int col1, int col2, int col3, double size)
{

    if(!alloc || update_needed)
	{
	    create_copy();

     	}

    if(col1 >= data->size2 || col2 >= data->size2 || col3 >= data->size2)
	{
	    throw invalid_argument("Invalid column selection for DS9 write");
	}

    ofstream output(filename.c_str());

    output <<"# Region file format: DS9 version 4.0" <<endl;
    output <<"# SaWLens catalogue writer: ASCII->DS9 conversion." <<endl;


    output <<"global color=white dashlist=8 3 width=1 font=\"helvetica 10 normal\" select=1 highlite=0 dash=0 fixed=0 edit=0 move=0 delete=1 include=1 source=0" <<endl;

    if(WCS)
	{
	    output <<"fk5" <<endl;
	}
    else
	{
	    output <<"linear" <<endl;
	}

    bool colour = false;

    for(int i = 0; i < data->size1; i++)
	{
	    output <<"circle(" <<gsl_matrix_get(data,i,col1) <<"," <<gsl_matrix_get(data,i,col2) <<"," <<size <<") #text={" <<gsl_matrix_get(data,i,col3) <<"}" <<" color=" <<flush;

	    if(colour)
		{
		    output <<"red" <<flush;
		    if(i != data->size1-1)
			{
			    if((int) gsl_matrix_get(data,i,col3) != (int) gsl_matrix_get(data,i+1,col3))
				{
				    colour = false;
				}
			}
		}
	    else
		{
		    output <<"green" <<flush;
		    if(i != data->size1-1)
			{
			    if((int) gsl_matrix_get(data,i,col3) != (int) gsl_matrix_get(data,i+1,col3))
				{
				    colour = true;
				}
			}
		}

	    output <<" tag={SL} tag={sys" <<(int) gsl_matrix_get(data,i,col3) <<"}" <<endl;
	}
    output.close();

} 


void CatWriter::add_tocol(int col, double value)
{

    if(col < 0 || col >= num_columns)
	{
	    throw invalid_argument("Invalid column selection for add_tocol");
	}

    if(!alloc || update_needed)
	{
	    create_copy();

     	}

    for(int i = 0; i < num_rows; i++)
	{
	    gsl_matrix_set(data,i,col,gsl_matrix_get(data,i,col)+value);
	}
}

void CatWriter::multiply_tocol(int col, double value)
{

    if(col < 0 || col >= num_columns)
	{
	    throw invalid_argument("Invalid column selection for multiply_tocol");
	}

    if(!alloc || update_needed)
	{
	    create_copy();

     	}

    for(int i = 0; i < num_rows; i++)
	{
	    gsl_matrix_set(data,i,col,gsl_matrix_get(data,i,col)*value);
	}
}

void CatWriter::power_tocol(int col, double value)
{

    if(col < 0 || col >= num_columns)
	{
	    throw invalid_argument("Invalid column selection for power_tocol");
	}

    if(!alloc || update_needed)
	{
	    create_copy();

     	}

    for(int i = 0; i < num_rows; i++)
	{
	    gsl_matrix_set(data,i,col,pow(gsl_matrix_get(data,i,col),value));
	}
}

void CatWriter::trig_tocol(int col, string selection)
{

    if(col < 0 || col >= num_columns)
	{
	    throw invalid_argument("Invalid column selection for trig_tocol");
	}

    if(!alloc || update_needed)
	{
	    create_copy();

     	}

    for(int i = 0; i < num_rows; i++)
	{
	    if(selection == "sin")
		{
		    gsl_matrix_set(data,i,col,sin(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "cos")
		{
		    gsl_matrix_set(data,i,col,cos(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "tan")
		{
		    gsl_matrix_set(data,i,col,tan(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "asin")
		{
		    gsl_matrix_set(data,i,col,asin(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "acos")
		{
		    gsl_matrix_set(data,i,col,acos(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "atan")
		{
		    gsl_matrix_set(data,i,col,atan(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "sinh")
		{
		    gsl_matrix_set(data,i,col,sinh(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "cosh")
		{
		    gsl_matrix_set(data,i,col,cosh(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "tanh")
		{
		    gsl_matrix_set(data,i,col,tanh(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "exp")
		{
		    gsl_matrix_set(data,i,col,exp(gsl_matrix_get(data,i,col)));
		}
	    else
		{
		    throw invalid_argument("Invalid selection for trig_tocol");
		}
	}
}

void CatWriter::addfunccol_tocol(int target, int col)
{

    if(target < 0 || col < 0 || target > num_columns || col > num_columns)
	{
	    throw invalid_argument("Invalid column selection for addfunc");
	}

    for(int i = 0; i < data->size1; i++)
	{
	    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)+gsl_matrix_get(data,i,col));
	}
}

void CatWriter::addfunccol_tocol(int target, int col, double power)
{

    if(target < 0 || col < 0 || target > num_columns || col > num_columns)
	{
	    throw invalid_argument("Invalid column selection for addfunc");
	}

    for(int i = 0; i < data->size1; i++)
	{
	    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)+pow(gsl_matrix_get(data,i,col),power));
	}
}

void CatWriter::addfunccol_tocol(int target, int col, string selection)
{

    if(target < 0 || col < 0 || target > num_columns || col > num_columns)
	{
	    throw invalid_argument("Invalid column selection for addfunc");
	}

    for(int i = 0; i < data->size1; i++)
	{
	    if(selection == "sin")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)+sin(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "cos")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)+cos(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "tan")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)+tan(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "asin")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)+asin(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "acos")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)+acos(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "atan")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)+atan(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "sinh")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)+sinh(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "cosh")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)+cosh(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "tanh")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)+tanh(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "exp")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)+tanh(gsl_matrix_get(data,i,col)));
		}
	    else
		{
		    throw invalid_argument("Invalid selection for addfunctocol");
		}
	}
}

void CatWriter::multiplyfunccol_tocol(int target, int col)
{

    if(target < 0 || col < 0 || target > num_columns || col > num_columns)
	{
	    throw invalid_argument("Invalid column selection for addfunc");
	}

    for(int i = 0; i < data->size1; i++)
	{
	    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)*gsl_matrix_get(data,i,col));
	}
}

void CatWriter::multiplyfunccol_tocol(int target, int col, double power)
{

    if(target < 0 || col < 0 || target > num_columns || col > num_columns)
	{
	    throw invalid_argument("Invalid column selection for addfunc");
	}

    for(int i = 0; i < data->size1; i++)
	{
	    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)*pow(gsl_matrix_get(data,i,col),power));
	}
}

void CatWriter::multiplyfunccol_tocol(int target, int col, string selection)
{

    if(target < 0 || col < 0 || target > num_columns || col > num_columns)
	{
	    throw invalid_argument("Invalid column selection for addfunc");
	}

    for(int i = 0; i < data->size1; i++)
	{
	    if(selection == "sin")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)*sin(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "cos")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)*cos(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "tan")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)*tan(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "asin")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)*asin(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "acos")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)*acos(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "atan")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)*atan(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "sinh")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)*sinh(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "cosh")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)*cosh(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "tanh")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)*tanh(gsl_matrix_get(data,i,col)));
		}
	    else if(selection == "exp")
		{
		    gsl_matrix_set(data,i,target,gsl_matrix_get(data,i,target)*tanh(gsl_matrix_get(data,i,col)));
		}
	    else
		{
		    throw invalid_argument("Invalid selection for addfunctocol");
		}
	}
}



void CatWriter::kill_row(int col, double thresh1, double thresh2)
{

    if(col < 0 || col >= num_columns)
	{
	    throw invalid_argument("Invalid column selection for trig_tocol");
	}

    if(!alloc || update_needed)
	{
	    create_copy();
     	}



    vector<bool> marker;
    int counter = 0;

    for(int i = 0; i < data->size1; i++)
	{
	    if(gsl_matrix_get(data,i,col) < thresh1 || gsl_matrix_get(data,i,col) > thresh2)
		{
		    marker.push_back(true);
		    counter++;
		}
	    else
		{
		    marker.push_back(false);
		}
	}
    gsl_matrix *nofuture = gsl_matrix_calloc(data->size1, data->size2);
    gsl_matrix_memcpy(nofuture,data);
    gsl_matrix_free(data);
    data = gsl_matrix_calloc(nofuture->size1-counter,nofuture->size2);

    int row_counter = 0;

    for(int i = 0; i < nofuture->size1; i++)
	{
	    if(!marker[i])
		{
		    for(int j = 0; j < nofuture->size2; j++)
			{
			    gsl_matrix_set(data,row_counter,j,gsl_matrix_get(nofuture,i,j));
			    row_counter++;
			}
		}
	}

    gsl_matrix_free(nofuture);
    num_rows = data->size1;
    marker.clear();
}


void CatWriter::add_col(int col, double value, string name)
{

    if(col < 0 || col > num_columns)
	{
	    throw invalid_argument("Invalid column selection for trig_tocol");
	}

    if(!alloc || update_needed)
	{
	    create_copy();
     	}

    gsl_matrix *nofuture = gsl_matrix_calloc(data->size1, data->size2);
    gsl_matrix_memcpy(nofuture,data);
    gsl_matrix_free(data);
    data = gsl_matrix_calloc(nofuture->size1,nofuture->size2+1);


    for(int i = 0; i < nofuture->size1; i++)
	{
	    bool gotta =  false;

	    for(int j = 0; j < nofuture->size2+1; j++)
		{
		    if(j == col)
			{
			    gsl_matrix_set(data,i,j,value);
			    gotta = true;
			}
		    else if(!gotta)
			{
			    gsl_matrix_set(data,i,j,gsl_matrix_get(nofuture,i,j));
			}
		    else
			{
			    gsl_matrix_set(data,i,j,gsl_matrix_get(nofuture,i,j-1));
			}

		}
	    
	}

    gsl_matrix_free(nofuture);
    num_columns = data->size2;

    vector<std::string> littlehelper;
    for(int i = 0; i < num_columns; i++)
	{
	    if(i < col)
		{
		    littlehelper.push_back(names[i]);
		}
	    else if(i == col)
		{
		    littlehelper.push_back(name);
		}
	    else
		{
		    littlehelper.push_back(names[i-1]);
		}
	}
    names.clear();
    names = littlehelper;
    littlehelper.clear();



}




















