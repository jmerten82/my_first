#include <iostream>
#include <ctime>
#include <string>
#include <sstream>
#include "options.h"
#include "mpi_inputfield.h"
#include "mpi.h"

using namespace std;

int main(int argc, char *argv[])

{
  int my_rank;
  int p;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  if(my_rank == 0)
    {

      cout <<"-------- SAW 1.4 Field preparation routine --------" <<endl;
      cout <<endl;
    }
      
  //Reading configuration file into options class
  FieldOptions options1(argv[1]);
  //taking the time

  long starttime;
  long stoptime;

  time(&starttime);


  //Starting the main loop through all resolutions

  for(int i = 0; i < options1.numres(); i++)
    {
      if(my_rank == 0)
	{
	  cout <<endl;
	  cout <<"x-dimension: " <<options1.showstep(i) <<endl;
	}

      DataGrid grid1(options1,i,my_rank,p);
      grid1.analyse(options1,my_rank,p);
      grid1.write(options1,i,my_rank);
      if(my_rank == 0)
	{
	  cout <<endl;
	}
    }

  time(&stoptime);

  if(my_rank == 0)
    {
      cout <<"Total runtime: " <<stoptime-starttime <<endl;
    }


  MPI_Finalize();
  return 0;

}

 

 



/*


int main(int argc, char *argv[])
{
  cout <<"----------SAW 1.0 Field preparation routine----------" <<endl;
  cout <<endl;;
  cout <<"Chosen configuration file: " <<argv[1] <<endl;

  //reading the config file and defining the options 

  FieldOptions options1(argv[1]);

  //declare some filenames and timimg variables


  string outfile1;
  string outfile2;
  string outfile3;
  string outfile4;
  string outfile5;
  string outfile6;
  string outfile7;
  string outfile8;
  long starttime;
  long stoptime;

  //starting the counter 
  
  time(&starttime);

  //going through the resolution loop

  for(int i = options1.startdim(); i <= options1.stopdim(); i++)
    {



      // defining the filenames for each resolution

      ostringstream output1;
      ostringstream output2;
      ostringstream output3;
      ostringstream output4;
      ostringstream output5;
      ostringstream output6;
      ostringstream output7;
      ostringstream output8;
      output1 <<options1.filename(2) <<i <<".dat";
      output2 <<options1.filename(2) <<i <<".fits";
      output3 <<options1.filename(5) <<i <<".dat";
      output4 <<options1.filename(5) <<i <<".fits";
      output5 <<options1.filename(6) <<i <<".dat";
      output6 <<options1.filename(6) <<i <<".fits";
      output7 <<options1.filename(9) <<i <<".dat";
      output8 <<options1.filename(9) <<i <<".fits";
      outfile1 = output1.str();
      outfile2 = output2.str();
      outfile3 = output3.str();
      outfile4 = output4.str();
      outfile5 = output5.str();
      outfile6 = output6.str();
      outfile7 = output7.str();
      outfile8 = output8.str();

      //constructing the major class from the filenames and options      

      
      WeightedInputField field1(options1.filename(3),options1.filename(0),i,options1.cutpoint(0),options1.cutpoint(1),options1.cutpoint(2),options1.cutpoint(3),options1.cutpoint(4),options1.cutpoint(5),options1.cutpoint(6),options1.cutpoint(7),options1.showflag(4),options1.galaxies(2));

      //masking the field, non tested for strong lensing regime
	       
      field1.checkmask(options1.filename(7));
      
      
      if(options1.showflag(0))
	{

	  //version for shear 
	  field1.average(options1.galaxies(0),options1.increment());
	  field1.docovariance(options1.covradius());
	  field1.showinfo(0);
	  field1.showinfo(0,outfile1);
	  field1.writetofits(outfile2,0);
	}
      if(options1.showflag(2))
	{

	  //version for the critical curves
	  field1.ccurvebuild(options1.filename(1));
	  field1.showinfo(1);
	  field1.showinfo(1,outfile5);
	  field1.writetofits(outfile6,1);
	}

      cout <<"Done!" <<endl;
      if(options1.showflag(3))
	{

	  //version for the multiple image systems
	  field1.msystembuild(options1.filename(8),outfile7);
	  field1.showinfo(1);
	}
      cout <<"Done!" <<endl;
      if(options1.showflag(1))
	{

	  //version for the flexion field, not very elegant yet
	  FlexionField flexion1(options1.filename(3),options1.filename(4),i,options1.galaxies(1),options1.cutpoint(0),options1.cutpoint(1),options1.cutpoint(2),options1.cutpoint(3));
	  
	  flexion1.average(options1.increment());
	  flexion1.docovariance_fast(options1.covradius());
	  flexion1.showinfo();
	  flexion1.showinfo(outfile3);
	  flexion1.writetofits(outfile4);
	}

      cout <<endl;
		
      outfile1 ="";
      outfile2 ="";
      outfile3 ="";
      outfile4 ="";
      outfile5 ="";
      outfile6 ="";
      outfile7 ="";
      outfile8 ="";
    }
 
  //stopping the stopwatch  
      
      time(&stoptime);
      
      cout <<"Total runtime: " <<stoptime-starttime <<endl;
    
*/
  

    
