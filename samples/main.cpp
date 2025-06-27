//UNN
//Litvyakov D. D.
//

#include "template.h"
#include <iostream>
#include "mpi.h"
#define type int

size_t msize = 2;


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

    type* lmatrix = new type[msize * msize];
    type* rmatrix = new type[msize * msize];
    type* resmatrix = new type[msize * msize];
	
    generatematrix(lmatrix, msize);
    generatematrix(rmatrix, msize);

    mul(lmatrix, rmatrix, resmatrix, msize);

    std::cout << checkmatrix(resmatrix,resmatrix,msize);





    delete [] lmatrix;
    delete [] rmatrix;
    delete [] resmatrix;
	MPI_Finalize();
	return 0;
}
