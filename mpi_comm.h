#ifndef   	MPI_COMM_H_
# define   	MPI_COMM_H_

#include <iostream>
#include <cmath> 
#include <gsl/gsl_matrix.h>
#include "mpi.h"

using namespace std;

/** 
    Routines to send and receive GSL arrays in an MPI environment.
**/

class MPI_Array_Manager
{

  /**

     class which manages the sizes and indexes in the matrices and vectors
     calculated on different processes and then send to a leading process.
     It makes sure that the right amount of memory is allocated and that 
     the order in the final matrix is correct.
     All the index-splitting is based on the row-index, columns in matrices
     are calculated on all processes.

  **/

 private:

  int my_rank;
  /*
    Rank of the process the class is executed.
  */
  int p;
  /* 
     Total number of processes.
  */

  bool creation;
  /*
    This process takes part in the matrix/vector-splitting process?
  */

  int size;
  /*
    The number of (row-)indices calculated on this process.
  */

  int totalsize;
  /*
    The total number of (row-)indices for the split matrix/vector.
  */

  int startindex;
  /*
    The (row-)index in the big matrix where this specific process
    array starts.
  */
  int stopindex;
  /*
    The (row-)index in the big matrix where the summation stops.
  */  

  bool root;
  /*
    Am I the root process which is meant to carry the full matrix in 
    the end?
  */

  int rootpos;
  /*
    Process number of the root process.
  */

  gsl_vector_int *array_processes;
  /*
    Vector containg the process ranks of the array processes.
  */

  int array_split_size;
  /*
    Number of prcoesses taking part in the splitting.
  */

  MPI_Comm array_comm;
  /*
    The MPI-communicator containing the processes in our matrix-splitting.
  */


 public:

  MPI_Array_Manager(int my_rankgiven, int pgiven, int x_dim, int y_dim, int rootgiven);
  /*
    Constructor, needs the MPI-Env. info, the x -and y-dimension if it is 
    a map which has to be splitted and the number of the root process.
  */

  MPI_Array_Manager(int my_rankgiven, int pgiven, int dim, int rootgiven);
  /*
    As above but with the absolute size of the index to be splitted.
  */
  ~MPI_Array_Manager();
  /*
    Standard destructor.
  */
  bool greenlight();
  /*
    Returns the flag if the process is part of the index-splitting.
  */
  int showsize();
  /*
    Returns the number of index elements on this process 0 if process is not
    part of the index splitting.
  */ 
  int showtotalsize();
  /*
    The total number of elements in the index which is split. 
  */ 
  bool showroot();
  /*
    Flag which is true if the process is the root process.
  */
  int showrootpos();
  /*
    Returns WORLD rank of the root process.
  */
  int showstartindex();
  /*
    Returns the original position of the first element of this process wrt
    the original index.
  */
  int showstopindex();
  /*
    Original position of the last element. 
  */
  void createcomm();
  /*
    Creates an MPI-communicator for the index-splitting task.
  */
  int showcommsize();
  /*
    Return the number of processes in the splitting communicator.
  */
  MPI_Comm arraycomm();
  /*
    Returns the splitting communicator.
  */
  void printarrayranks();
  /*
    Prints the processes, taking part in the splitting to stdout.
  */  

};

void recv_gsl_fromworld_new(gsl_matrix *inmatrix, gsl_matrix *outmatrix);
/*
  Sends an gsl_matrix to MPI_COMM and adds all matrices at the process
  with my_rank 0. Equivalent to a MPI_Reduce.
*/

void recv_gsl_fromworld_new(gsl_matrix_int *inmatrix, gsl_matrix_int *outmatrix);
/*
  Sends an gsl_matrix to MPI_COMM and adds all matrices at the process
  with my_rank 0. Equivalent to a MPI_Reduce.
*/
void recv_gsl_fromworld_new(gsl_vector *invector, gsl_vector *outvector);
/*
  Sends an gsl_vector to MPI_COMM and adds all matrices at the process
  with my_rank 0. Equivalent to a MPI_Reduce.
*/
void recv_gsl_fromworld_new(gsl_vector_int *invector, gsl_vector_int *outvector);
/*
  Sends an gsl_vector to MPI_COMM and adds all matrices at the process
  with my_rank 0. Equivalent to a MPI_Reduce.
*/

void recv_gsl_fromworld_new(gsl_matrix *outmatrix);
/*
  As above but allows aliasing
*/

void recv_gsl_fromworld_new(gsl_matrix_int *outmatrix);

void recv_gsl_fromworld_new(gsl_vector *outvector);

void recv_gsl_fromworld_new(gsl_vector_int *outvector);

void recv_gsl_fromworld_ultranew(gsl_matrix *outmatrix,int my_rank,int p);

void recv_gsl_fromworld_ultranew(gsl_matrix_int *outmatrix,int my_rank,int p);

void recv_gsl_fromworld_ultranew(gsl_vector *outvector,int my_rank,int p);

void recv_gsl_fromworld_ultranew(gsl_vector_int *outvector,int my_rank,int p);



void send_gsl_toworld(gsl_matrix *matrix);
/*
  Sends a gsl_matrix from process with my_rank 0 to all other processes.
  Equivalent to MPI_Bcast.
*/

void send_gsl_toworld(gsl_matrix_int *matrix);
/*
  Sends a gsl_matrix from process with my_rank 0 to all other processes.
  Equivalent to MPI_Bcast.
*/

void send_gsl_toworld(gsl_matrix_uchar *matrix);
/*
  Sends a gsl_matrix from process with my_rank 0 to all other processes.
  Equivalent to MPI_Bcast.
*/

void send_gsl_toworld(gsl_vector *vector);
/*
  Sends a gsl_vector from process with my_rank 0 to all other processes.
  Equivalent to MPI_Bcast.
*/
void send_gsl_toworld(gsl_vector_int *vector);
/*
  Sends a gsl_vector from process with my_rank 0 to all other processes.
  Equivalent to MPI_Bcast.
*/
void send_gsl_toworld(gsl_vector_uchar *vector);
/*
  Sends a gsl_vector from process with my_rank 0 to all other processes.
  Equivalent to MPI_Bcast.
*/

void send_gsl_toworld(gsl_matrix *matrix, MPI_Array_Manager &array1);
/*
  As above, but broadcasting is just done within a certain Communicator,
  defined by the Array_Manager.
*/

void send_gsl_toworld(gsl_matrix_int *matrix, MPI_Array_Manager &array1);
void send_gsl_toworld(gsl_matrix_uchar *matrix, MPI_Array_Manager &array1);
void send_gsl_toworld(gsl_vector *vector, MPI_Array_Manager &array1);
void send_gsl_toworld(gsl_vector_int *vector, MPI_Array_Manager &array1);
void send_gsl_toworld(gsl_vector_uchar *vector, MPI_Array_Manager &array1);

void recv_gsl_fromworld(gsl_matrix *sendmatrix, gsl_matrix *recvmatrix,MPI_Array_Manager &array1);

/*
  Receives the data in sendmatrix and stores in rank-order in recvmatrix on 
  process root. Process rank and total number of processes have to be given.
*/

void recv_gsl_fromworld(gsl_matrix_int *sendmatrix, gsl_matrix_int *recvmatrix,MPI_Array_Manager &array1);
void recv_gsl_fromworld(gsl_matrix_uchar *sendmatrix, gsl_matrix_uchar *recvmatrix,MPI_Array_Manager &array1);

void recv_gsl_fromworld(gsl_vector *sendvector, gsl_vector *recvvector,MPI_Array_Manager &array1);

void recv_gsl_fromworld(gsl_vector_int *sendvector, gsl_vector_int *recvvector,MPI_Array_Manager &array1);
void recv_gsl_fromworld(gsl_vector_uchar *sendvector, gsl_vector_uchar *recvvector,MPI_Array_Manager &array1);

void all_recv_gsl_fromworld(gsl_matrix *sendmatrix, gsl_matrix *recvmatrix,MPI_Array_Manager &array1);
/*
  Refers to Allgather, not working properly.
*/

void all_recv_gsl_fromworld(gsl_matrix_int *sendmatrix, gsl_matrix_int *recvmatrix,MPI_Array_Manager &array1);
void all_recv_gsl_fromworld(gsl_matrix_uchar *sendmatrix, gsl_matrix_uchar *recvmatrix,MPI_Array_Manager &array1);

void all_recv_gsl_fromworld(gsl_vector *sendvector, gsl_vector *recvvector,MPI_Array_Manager &array1);

void all_recv_gsl_fromworld(gsl_vector_int *sendvector, gsl_vector_int *recvvector,MPI_Array_Manager &array1);
void all_recv_gsl_fromworld(gsl_vector_uchar *sendvector, gsl_vector_uchar *recvvector,MPI_Array_Manager &array1);

void special_recv_gsl_fromworld(gsl_matrix *sendmatrix,gsl_matrix *recvmatrix, MPI_Array_Manager &array1);
/*
  Allgather within an special communicator defined by the Array Manager. Is
  not working properly.
*/

void special_recv_gsl_fromworld(gsl_matrix_int *sendmatrix,gsl_matrix_int *recvmatrix, MPI_Array_Manager &array1);
void special_recv_gsl_fromworld(gsl_matrix_uchar *sendmatrix,gsl_matrix_uchar *recvmatrix, MPI_Array_Manager &array1);

void special_recv_gsl_fromworld(gsl_vector *sendmatrix,gsl_vector *recvmatrix, MPI_Array_Manager &array1);

void special_recv_gsl_fromworld(gsl_vector_int *sendmatrix,gsl_vector_int *recvmatrix, MPI_Array_Manager &array1);
void special_recv_gsl_fromworld(gsl_vector_uchar *sendmatrix,gsl_vector_uchar *recvmatrix, MPI_Array_Manager &array1);


#endif 	    /* !MPI_COMM_H_ */
