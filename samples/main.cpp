//UNN
//Litvyakov D. D.
//

#include "template.h"
#include <iostream>
#include <chrono>
#define type int

size_t msize = 1024;


int main(int argc, char** argv)
{
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();

    int myid, numprocs;
    MPI_Status status;
    int  namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Get_processor_name(processor_name, &namelen);

    type* lmatrix = nullptr;
    type* rmatrix = nullptr;
    type* resmatrix1 = nullptr;
    type* resmatrix2 = nullptr;
    type* resmatrix3 = nullptr;

    if (myid == 0) {

        lmatrix = new type[msize * msize];
        rmatrix = new type[msize * msize];
        resmatrix1 = new type[msize * msize];
        resmatrix2 = new type[msize * msize];
        resmatrix3 = new type[msize * msize];

        generatematrix(lmatrix, msize);
        generatematrix(rmatrix, msize);

        start = std::chrono::steady_clock::now();
        basedmul(lmatrix, rmatrix, resmatrix1, msize);
        end = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
        std::cout << "standart mul time: " << elapsed / 1000. << " seconds" << std::endl;

        start = std::chrono::steady_clock::now();
        mul(lmatrix, rmatrix, resmatrix2, msize);
        end = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
        std::cout << "cache friendly mul time: " << elapsed / 1000. << " seconds" << std::endl;


        std::cout << "Result is correct? " << checkmatrix(resmatrix1, resmatrix2, msize) << std::endl;

        start = std::chrono::steady_clock::now();
    }

    MatrixMultiplicationMPI(lmatrix, rmatrix, resmatrix3, msize);

        if (myid == 0) {
            end = std::chrono::steady_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
            std::cout << "mpi time: " << elapsed / 1000. << " seconds" << std::endl;


            std::cout << "Result is correct? " << checkmatrix(resmatrix1, resmatrix3, msize) << std::endl;

            delete[] lmatrix;
            delete[] rmatrix;
            delete[] resmatrix1;
            delete[] resmatrix2;
            delete[] resmatrix3;
            lmatrix = rmatrix = resmatrix1 = resmatrix2 = resmatrix3 = nullptr;
        }

    MPI_Finalize();

	return 0;
}
