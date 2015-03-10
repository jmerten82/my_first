#ifndef   	MPI_COREROUTINES_H_
# define   	MPI_COREROUTINES_H_

#include <iostream>
#include <cmath>
#include <gsl/gsl_matrix.h>
#include "options.h"
#include "galaxycluster.h"
#include "util.h"
#include "mpi.h"
#include "rw_fits.h"
#include "mpi_comm.h"
#include "mpi_datacontainer.h"

using namespace std;



void MPIWeakFactor(ReconstructionOptions &options,GalaxyCluster &referencecluster, MPIDataContainer &container);
/*
  Calcluates the prefactor matrices for shear and flexion which are Fij*Zi*Zj
  in the notation of arXiv:0806.1967.
*/

void MPIStrongFactor(ReconstructionOptions &options,GalaxyCluster &referencecluster, MPIDataContainer &container);
/*
  Calculates the prefactor matrix for the critical curve term as
  in arXiv:0806.1967.
*/


void MPI_Blk(ReconstructionOptions &options,GalaxyCluster &referencecluster,MPIDataContainer &container,int p, int my_rank);
/*
  Constructs the big reconstruction coefficient matrix. Depending on the 
  ReconstructionOptions this matrix takes into account shear, flexion, critical
  curve and multiple image systems. The referencecluster delivers the
  regularisation values.
  Depending on the rank of the process and the total number of processes, 
  the matrix building is split between the processes. 
*/

void MPI_Vl(ReconstructionOptions &options,GalaxyCluster &referencecluster,MPIDataContainer &container,int p, int my_rank);
/*
  Constructs the big reconstruction data vector. Depending on the 
  ReconstructionOptions this matrix takes into account shear, flexion, critical
  curve and multiple image systems. The referencecluster delivers the
  regularisation values. 
  Depending on the rank of the process and the total number of processes, 
  the matrix building is split between the processes. 
*/



#endif 	    /* !MPI_COREROUTINES_H_ */
