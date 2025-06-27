//UNN
//Litvyakov D. D.
//

#include "template.h"
#include <iostream>
#include "mpi.h"

int main(int argc, char** argv)
{
    int myid, numprocs;
    MPI_Status status;
    int  namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Get_processor_name(processor_name, &namelen);


	std::cout<<check();


	MPI_Finalize();
	return 0;
}
