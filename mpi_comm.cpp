#include "mpi_comm.h"


MPI_Array_Manager::MPI_Array_Manager(int my_rankgiven, int pgiven, int x_dim, int y_dim, int rootgiven)
{

  my_rank = my_rankgiven;
  p = pgiven;

  totalsize = x_dim*y_dim;

  int frac = (int) ceil((double) totalsize/p);

  if((my_rank+1)*frac <= totalsize)
    {
      creation = true;
      size = frac;
      startindex = my_rank*frac;
      stopindex = startindex+size-1;
    }
  else if(my_rank*frac < totalsize)
    {
      creation = true;
      size = totalsize - my_rank*frac;
      startindex = my_rank*frac;
      stopindex = startindex+size-1;
    }
  else
    {
      size = 0;
      creation = false;
      startindex = -1;
      stopindex = -2;
    }

  if(my_rank == rootgiven)
    {
      root = true;
    }
  else
    {
      root = false;
    }

  rootpos = rootgiven;

  int counter = 0;
  for(int i = 0; i < p; i++)
    {
      if(i*frac < totalsize)
	{
	  counter++;
	}
    }
  array_processes = gsl_vector_int_calloc(counter);
  array_split_size = counter;
  for(int i = 0; i < array_split_size;i++)
    {
      gsl_vector_int_set(array_processes,i,i);
    }
  
}
MPI_Array_Manager::MPI_Array_Manager(int my_rankgiven, int pgiven, int dim, int rootgiven)
{

  my_rank = my_rankgiven;
  p = pgiven;
  totalsize = dim;
  int frac = (int) ceil((double) totalsize/p);

  if((my_rank+1)*frac <= totalsize)
    {
      creation = true;
      size = frac;
      startindex = my_rank*frac;
      stopindex = startindex+size-1;
    }
  else if(my_rank*frac < totalsize)
    {
      creation = true;
      size = totalsize - my_rank*frac;
      startindex = my_rank*frac;
      stopindex = startindex+size-1;
    }
  else
    {
      size = 0;
      creation = false;
      startindex = -1;
      stopindex = -2;
    }

  if(my_rank == rootgiven)
    {
      root = true;
    }
  else
    {
      root = false;
    }

 rootpos = rootgiven;
 int counter = 0;
 for(int i = 0; i < p; i++)
   {
     if(i*frac < totalsize)
       {
	 counter++;
       }
   }
 array_processes = gsl_vector_int_calloc(counter);
 array_split_size = counter;
 for(int i = 0; i < array_split_size;i++)
   {
     gsl_vector_int_set(array_processes,i,i);
   }
	  
}

MPI_Array_Manager::~MPI_Array_Manager()
{
  creation = false;
  root = 0;
}

bool MPI_Array_Manager::greenlight()
{
  return creation;
}

int MPI_Array_Manager::showsize()
{

  return size;
}

int MPI_Array_Manager::showtotalsize()
{
  return totalsize;
}

bool MPI_Array_Manager::showroot()
{
  return root;
}
int MPI_Array_Manager::showrootpos()
{
  return rootpos;
}

int MPI_Array_Manager::showstartindex()
{
  return startindex;
}

int MPI_Array_Manager::showstopindex()
{
  return stopindex;
}

void MPI_Array_Manager::createcomm()
{
  MPI_Group MPI_GROUP_WORLD;
  MPI_Group array_group;

  MPI_Comm_group(MPI_COMM_WORLD,&MPI_GROUP_WORLD);

  MPI_Group_incl(MPI_GROUP_WORLD,array_split_size,gsl_vector_int_ptr(array_processes,0),&array_group);

  MPI_Comm_create(MPI_COMM_WORLD,array_group,&array_comm);

}

MPI_Comm MPI_Array_Manager::arraycomm()
{
  return array_comm;
}

int MPI_Array_Manager::showcommsize()
{
  return array_split_size;
}

void MPI_Array_Manager::printarrayranks()
{
  for(int i = 0; i < array_processes->size; i++)
    {
      cout <<gsl_vector_int_get(array_processes,i) <<"\t" <<flush;
    }
  cout <<endl;
}


void send_gsl_toworld(gsl_matrix *matrix)
{
  MPI_Bcast(gsl_matrix_ptr(matrix,0,0),matrix->size1*matrix->size2,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void send_gsl_toworld(gsl_matrix_int *matrix)
{
  MPI_Bcast(gsl_matrix_int_ptr(matrix,0,0),matrix->size1*matrix->size2,MPI_INT,0,MPI_COMM_WORLD);
}

void send_gsl_toworld(gsl_matrix_uchar *matrix)
{
  MPI_Bcast(gsl_matrix_uchar_ptr(matrix,0,0),matrix->size1*matrix->size2,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);
}

void send_gsl_toworld(gsl_vector *vector)
{
  MPI_Bcast(gsl_vector_ptr(vector,0),vector->size,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void send_gsl_toworld(gsl_vector_int *vector)
{
  MPI_Bcast(gsl_vector_int_ptr(vector,0),vector->size,MPI_INT,0,MPI_COMM_WORLD);
}

void send_gsl_toworld(gsl_vector_uchar *vector)
{
  MPI_Bcast(gsl_vector_uchar_ptr(vector,0),vector->size,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);
}

void send_gsl_toworld(gsl_matrix *matrix, MPI_Array_Manager &array1)
{
  MPI_Bcast(gsl_matrix_ptr(matrix,0,0),matrix->size1*matrix->size2,MPI_DOUBLE,0,array1.arraycomm());
}
void send_gsl_toworld(gsl_matrix_int *matrix, MPI_Array_Manager &array1)
{
  MPI_Bcast(gsl_matrix_int_ptr(matrix,0,0),matrix->size1*matrix->size2,MPI_INT,0,array1.arraycomm());
}
void send_gsl_toworld(gsl_matrix_uchar *matrix, MPI_Array_Manager &array1)
{
  MPI_Bcast(gsl_matrix_uchar_ptr(matrix,0,0),matrix->size1*matrix->size2,MPI_UNSIGNED_CHAR,0,array1.arraycomm());
}
void send_gsl_toworld(gsl_vector *vector,MPI_Array_Manager &array1)
{
  MPI_Bcast(gsl_vector_ptr(vector,0),vector->size,MPI_DOUBLE,0,array1.arraycomm());
}
void send_gsl_toworld(gsl_vector_int *vector,MPI_Array_Manager &array1)
{
  MPI_Bcast(gsl_vector_int_ptr(vector,0),vector->size,MPI_INT,0,array1.arraycomm());
}
void send_gsl_toworld(gsl_vector_uchar *vector,MPI_Array_Manager &array1)
{
  MPI_Bcast(gsl_vector_uchar_ptr(vector,0),vector->size,MPI_UNSIGNED_CHAR,0,array1.arraycomm());
}

void recv_gsl_fromworld_new(gsl_matrix *inmatrix, gsl_matrix *outmatrix)
{
  MPI_Reduce(inmatrix->data,outmatrix->data, inmatrix->size1*inmatrix->size2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
}
void recv_gsl_fromworld_new(gsl_matrix_int *inmatrix, gsl_matrix_int *outmatrix)
{
  MPI_Reduce(inmatrix->data,outmatrix->data, inmatrix->size1*inmatrix->size2,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
}
void recv_gsl_fromworld_new(gsl_vector *invector, gsl_vector *outvector)
{
  MPI_Reduce(invector->data,outvector->data, invector->size,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
}
void recv_gsl_fromworld_new(gsl_vector_int *invector, gsl_vector_int *outvector)
{
  MPI_Reduce(invector->data,outvector->data, invector->size,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
}
void recv_gsl_fromworld_new(gsl_matrix *outmatrix)
{
  gsl_matrix *nofuture = gsl_matrix_calloc(outmatrix->size1,outmatrix->size2);
  gsl_matrix_memcpy(nofuture,outmatrix);
  MPI_Reduce(nofuture->data,outmatrix->data,outmatrix->size1*outmatrix->size2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  gsl_matrix_free(nofuture);
}

void recv_gsl_fromworld_new(gsl_matrix_int *outmatrix)
{
  gsl_matrix_int *nofuture = gsl_matrix_int_calloc(outmatrix->size1,outmatrix->size2);
  gsl_matrix_int_memcpy(nofuture,outmatrix);
  MPI_Reduce(nofuture->data,outmatrix->data, outmatrix->size1*outmatrix->size2,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  gsl_matrix_int_free(nofuture);
}

void recv_gsl_fromworld_new(gsl_vector *outvector)
{
  gsl_vector *nofuture = gsl_vector_calloc(outvector->size);
  gsl_vector_memcpy(nofuture,outvector);
  MPI_Reduce(nofuture->data,outvector->data,outvector->size, MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  gsl_vector_free(nofuture);
}

void recv_gsl_fromworld_new(gsl_vector_int *outvector)
{
  gsl_vector_int *nofuture = gsl_vector_int_calloc(outvector->size);
  gsl_vector_int_memcpy(nofuture,outvector);
  MPI_Reduce(nofuture->data,outvector->data,outvector->size, MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  gsl_vector_int_free(nofuture);
}

void recv_gsl_fromworld_ultranew(gsl_matrix *outmatrix,int my_rank, int p)
{

  int frac = (int) ceil( (double) outmatrix->size1 / p);
  if((my_rank+1)*frac <= outmatrix->size1)
    {
      MPI_Gather(gsl_matrix_ptr(outmatrix,my_rank*frac,0),frac*outmatrix->size2, MPI_DOUBLE, gsl_matrix_ptr(outmatrix,0,0),frac*outmatrix->size2,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
  else if(my_rank*frac < outmatrix->size1)
    {
     MPI_Gather(gsl_matrix_ptr(outmatrix,my_rank*frac,0),(outmatrix->size1 - my_rank*frac)*outmatrix->size2, MPI_DOUBLE, gsl_matrix_ptr(outmatrix,0,0),(outmatrix->size1 - my_rank*frac)*outmatrix->size2, MPI_DOUBLE,0,MPI_COMM_WORLD);
    }

}

void recv_gsl_fromworld_ultranew(gsl_matrix_int *outmatrix,int my_rank, int p)
{

  int frac = (int) ceil( (double) outmatrix->size1 / p);
  if((my_rank+1)*frac <= outmatrix->size1)
    {
      MPI_Gather(gsl_matrix_int_ptr(outmatrix,my_rank*frac,0),frac*outmatrix->size2, MPI_INT, gsl_matrix_int_ptr(outmatrix,0,0),frac*outmatrix->size2,MPI_INT,0,MPI_COMM_WORLD);
    }
  else if(my_rank*frac < outmatrix->size1)
    {
     MPI_Gather(gsl_matrix_int_ptr(outmatrix,my_rank*frac,0),(outmatrix->size1 - my_rank*frac)*outmatrix->size2, MPI_INT, gsl_matrix_int_ptr(outmatrix,0,0),(outmatrix->size1 - my_rank*frac)*outmatrix->size2, MPI_INT,0,MPI_COMM_WORLD);
    }

}

void recv_gsl_fromworld_ultranew(gsl_vector *outvector,int my_rank, int p)
{

  int frac = (int) ceil( (double) outvector->size / p);
  if((my_rank+1)*frac <= outvector->size)
    {
      MPI_Gather(gsl_vector_ptr(outvector,my_rank*frac),frac, MPI_DOUBLE, gsl_vector_ptr(outvector,0),frac,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
  else if(my_rank*frac < outvector->size)
    {
     MPI_Gather(gsl_vector_ptr(outvector,my_rank*frac),(outvector->size - my_rank*frac), MPI_DOUBLE, gsl_vector_ptr(outvector,0),(outvector->size - my_rank*frac), MPI_DOUBLE,0,MPI_COMM_WORLD);
    }

}

void recv_gsl_fromworld_ultranew(gsl_vector_int *outvector,int my_rank, int p)
{

  int frac = (int) ceil( (double) outvector->size / p);
  if((my_rank+1)*frac <= outvector->size)
    {
      MPI_Gather(gsl_vector_int_ptr(outvector,my_rank*frac),frac, MPI_INT, gsl_vector_int_ptr(outvector,0),frac,MPI_INT,0,MPI_COMM_WORLD);
    }
  else if(my_rank*frac < outvector->size)
    {
     MPI_Gather(gsl_vector_int_ptr(outvector,my_rank*frac),(outvector->size - my_rank*frac), MPI_INT, gsl_vector_int_ptr(outvector,0),(outvector->size - my_rank*frac), MPI_INT,0,MPI_COMM_WORLD);
    }

}

void recv_gsl_fromworld(gsl_matrix *sendmatrix,gsl_matrix *recvmatrix, MPI_Array_Manager &array1)
{
  MPI_Gather(gsl_matrix_ptr(sendmatrix,0,0),sendmatrix->size1*sendmatrix->size2,MPI_DOUBLE,gsl_matrix_ptr(recvmatrix,0,0),sendmatrix->size1*sendmatrix->size2,MPI_DOUBLE,array1.showrootpos(),array1.arraycomm());
}

void recv_gsl_fromworld(gsl_matrix_int *sendmatrix,gsl_matrix_int *recvmatrix,MPI_Array_Manager &array1)
{
  MPI_Gather(gsl_matrix_int_ptr(sendmatrix,0,0),sendmatrix->size2*sendmatrix->size1,MPI_INT,gsl_matrix_int_ptr(recvmatrix,0,0),sendmatrix->size2*sendmatrix->size1,MPI_INT,array1.showrootpos(),array1.arraycomm());
}

void recv_gsl_fromworld(gsl_matrix_uchar *sendmatrix,gsl_matrix_uchar *recvmatrix,MPI_Array_Manager &array1)
{
  MPI_Gather(gsl_matrix_uchar_ptr(sendmatrix,0,0),sendmatrix->size2*sendmatrix->size1,MPI_UNSIGNED_CHAR,gsl_matrix_uchar_ptr(recvmatrix,0,0),sendmatrix->size2*sendmatrix->size1,MPI_UNSIGNED_CHAR,array1.showrootpos(),array1.arraycomm());
} 

void recv_gsl_fromworld(gsl_vector *sendvector,gsl_vector *recvvector,MPI_Array_Manager &array1)
{
  MPI_Gather(gsl_vector_ptr(sendvector,0),sendvector->size,MPI_DOUBLE,gsl_vector_ptr(recvvector,0),sendvector->size,MPI_DOUBLE,array1.showrootpos(),array1.arraycomm());
}

void recv_gsl_fromworld(gsl_vector_int *sendvector,gsl_vector_int *recvvector,MPI_Array_Manager &array1)
{
  MPI_Gather(gsl_vector_int_ptr(sendvector,0),sendvector->size,MPI_INT,gsl_vector_int_ptr(recvvector,0),sendvector->size,MPI_INT,array1.showrootpos(),array1.arraycomm());
}
void recv_gsl_fromworld(gsl_vector_uchar *sendvector,gsl_vector_uchar *recvvector,MPI_Array_Manager &array1)
{
  MPI_Gather(gsl_vector_uchar_ptr(sendvector,0),sendvector->size,MPI_UNSIGNED_CHAR,gsl_vector_uchar_ptr(recvvector,0),sendvector->size,MPI_UNSIGNED_CHAR,array1.showrootpos(),array1.arraycomm());
}

void special_recv_gsl_fromworld(gsl_matrix *sendmatrix,gsl_matrix *recvmatrix, MPI_Array_Manager &array1)
{
  if(array1.greenlight())
    {
      MPI_Allgather(gsl_matrix_ptr(sendmatrix,0,0),sendmatrix->size1*sendmatrix->size2,MPI_DOUBLE,gsl_matrix_ptr(recvmatrix,0,0),sendmatrix->size1*sendmatrix->size2,MPI_DOUBLE,MPI_COMM_WORLD);
    }
  else
    {
      MPI_Allgather(gsl_matrix_ptr(sendmatrix,0,0),0,MPI_DOUBLE,gsl_matrix_ptr(recvmatrix,0,0),0,MPI_DOUBLE,MPI_COMM_WORLD);
    }
}

void special_recv_gsl_fromworld(gsl_matrix_int *sendmatrix,gsl_matrix_int *recvmatrix, MPI_Array_Manager &array1)
{
  if(array1.greenlight())
    {
      MPI_Allgather(gsl_matrix_int_ptr(sendmatrix,0,0),sendmatrix->size1*sendmatrix->size2,MPI_INT,gsl_matrix_int_ptr(recvmatrix,0,0),sendmatrix->size1*sendmatrix->size2,MPI_INT,MPI_COMM_WORLD);
    }
  else
    {
      MPI_Allgather(gsl_matrix_int_ptr(sendmatrix,0,0),0,MPI_INT,gsl_matrix_int_ptr(recvmatrix,0,0),0,MPI_INT,MPI_COMM_WORLD);
    }
}

void special_recv_gsl_fromworld(gsl_vector *sendvector,gsl_vector *recvvector,MPI_Array_Manager &array1)
{

  if(array1.greenlight())
    {
      MPI_Allgather(gsl_vector_ptr(sendvector,0),sendvector->size,MPI_DOUBLE,gsl_vector_ptr(recvvector,0),sendvector->size,MPI_DOUBLE,MPI_COMM_WORLD);
    }
  else
    {
      MPI_Allgather(gsl_vector_ptr(sendvector,0),0,MPI_DOUBLE,gsl_vector_ptr(recvvector,0),0,MPI_DOUBLE,MPI_COMM_WORLD);
    }

}

void special_recv_gsl_fromworld(gsl_vector_int *sendvector,gsl_vector_int *recvvector,MPI_Array_Manager &array1)
{

  if(array1.greenlight())
    {
      MPI_Allgather(gsl_vector_int_ptr(sendvector,0),sendvector->size,MPI_INT,gsl_vector_int_ptr(recvvector,0),sendvector->size,MPI_INT,MPI_COMM_WORLD);
    }
  else
    {
      MPI_Allgather(gsl_vector_int_ptr(sendvector,0),0,MPI_INT,gsl_vector_int_ptr(recvvector,0),0,MPI_INT,MPI_COMM_WORLD);
    }

}


void all_recv_gsl_fromworld(gsl_matrix *sendmatrix,gsl_matrix *recvmatrix, MPI_Array_Manager &array1)
{
  MPI_Allgather(gsl_matrix_ptr(sendmatrix,0,0),sendmatrix->size1*sendmatrix->size2,MPI_DOUBLE,gsl_matrix_ptr(recvmatrix,0,0),sendmatrix->size1*sendmatrix->size2,MPI_DOUBLE,array1.arraycomm());
}

void all_recv_gsl_fromworld(gsl_matrix_int *sendmatrix,gsl_matrix_int *recvmatrix,MPI_Array_Manager &array1)
{
  MPI_Allgather(gsl_matrix_int_ptr(sendmatrix,0,0),sendmatrix->size2*sendmatrix->size1,MPI_INT,gsl_matrix_int_ptr(recvmatrix,0,0),sendmatrix->size2*sendmatrix->size1,MPI_INT,array1.arraycomm());
} 

void all_recv_gsl_fromworld(gsl_vector *sendvector,gsl_vector *recvvector,MPI_Array_Manager &array1)
{
  MPI_Allgather(gsl_vector_ptr(sendvector,0),sendvector->size,MPI_DOUBLE,gsl_vector_ptr(recvvector,0),sendvector->size,MPI_DOUBLE,array1.arraycomm());
}

void all_recv_gsl_fromworld(gsl_vector_int *sendvector,gsl_vector_int *recvvector,MPI_Array_Manager &array1)
{
  MPI_Allgather(gsl_vector_int_ptr(sendvector,0),sendvector->size,MPI_INT,gsl_vector_int_ptr(recvvector,0),sendvector->size,MPI_INT,array1.arraycomm());
}





