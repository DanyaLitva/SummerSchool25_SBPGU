//UNN
//Litvyakov D. D.
//

#pragma once
#include <random>
#include "mpi.h"
#include <iostream>
#include <vector>
#include <stdexcept>


int MinVal = -1000;
int MaxVal = 1000;


template <typename type>
void basedmul(type* left, type* right, type* result, size_t sz) {
    for (size_t i = 0; i < sz * sz; ++i) result[i] = type(0);

    for (size_t i = 0; i < sz; ++i) {
        for (size_t j = 0; j < sz; ++j) {
            for (size_t k = 0; k < sz; ++k) {
                result[i * sz + j] += left[i * sz + k] * right[k * sz + j];
            }
        }
    }
}


template <typename type>
void mul(type* left, type* right, type* result, size_t sz) {
    for (size_t i = 0; i < sz * sz; ++i) result[i] = type(0);

    for (size_t i = 0; i < sz; ++i) {
        type* lefttemp = &left[i * sz];
        type* resulttemp = &result[i * sz];
        for (size_t k = 0; k < sz; ++k) {  
            type tempval = lefttemp[k];            
            type* righttemp = &right[k * sz];
            for (size_t j = 0; j < sz; ++j) {
                resulttemp[j] += tempval *  righttemp[j];
            }
        }
    }
}


template <typename type>
bool checkmatrix(type* left, type* right, size_t sz) {
    for (size_t i = 0; i < (sz * sz); ++i) if (left[i] != right[i]) return false;
    return true;
}

template <typename type>
void generatematrix(type* matrix, size_t sz) {
    std::random_device r;
    std::default_random_engine e(r());
    std::uniform_real_distribution<double> coef_gen(MinVal,MaxVal);

    for (size_t i = 0; i < sz*sz; ++i) {
        matrix[i] = type(coef_gen(e));
    }
}






template <typename type>
void Flip(type* Matrix, size_t Size) {
    type temp;
    for (size_t i = 0; i < Size; ++i) {
        for (size_t j = i + 1; j < Size; ++j) {
            temp = Matrix[i * Size + j];
            Matrix[i * Size + j] = Matrix[j * Size + i];
            Matrix[j * Size + i] = temp;
        }
    }
}




//the algorithm from the unn website http://www.hpcc.unn.ru/?dir=910
template <typename type>
void MatrixMultiplicationMPI(type*& A, type*& B, type*& C, size_t& Size) {

    int ProcRank, ProcNum;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);

    if (Size % ProcNum != 0) {
        if (ProcRank == 0) {
            std::cerr << "Error: Matrix size must be divisible by number of processes\n";
        }
        MPI_Finalize();
        exit(1);
    }

    MPI_Datatype mpi_type;
    if (typeid(type) == typeid(double)) mpi_type = MPI_DOUBLE;
    else if (typeid(type) == typeid(float)) mpi_type = MPI_FLOAT;
    else if (typeid(type) == typeid(int)) mpi_type = MPI_INT;
    else throw std::runtime_error("Unsupported matrix type");

    int dim = Size;
    int i, j, k, p, ind;
    type temp;
    MPI_Status Status;
    int ProcPartSize = dim / ProcNum;
    int ProcPartElem = ProcPartSize * dim;
    type* bufA = new type[ProcPartElem];
    type* bufB = new type[ProcPartElem];
    type* bufC = new type[ProcPartElem];
    int ProcPart = dim / ProcNum, part = ProcPart * dim;
    if (ProcRank == 0) {
        Flip(B, Size);
    }
    


    MPI_Scatter(A, part, mpi_type, bufA, part, mpi_type, 0, MPI_COMM_WORLD);
    MPI_Scatter(B, part, mpi_type, bufB, part, mpi_type, 0, MPI_COMM_WORLD);

    temp = 0.0;
    for (i = 0; i < ProcPartSize; i++) {
        for (j = 0; j < ProcPartSize; j++) {
            for (k = 0; k < dim; k++)
                temp += bufA[i * dim + k] * bufB[j * dim + k];
            bufC[i * dim + j + ProcPartSize * ProcRank] = temp;
            temp = 0.0;
        }
    }

    int NextProc; int PrevProc;
    for (p = 1; p < ProcNum; p++) {
        NextProc = ProcRank + 1;
        if (ProcRank == ProcNum - 1)
            NextProc = 0;
        PrevProc = ProcRank - 1;
        if (ProcRank == 0)
            PrevProc = ProcNum - 1;
        MPI_Sendrecv_replace(bufB, part, mpi_type, NextProc, 0, PrevProc, 0, MPI_COMM_WORLD, &Status);
        temp = 0.0;
        for (i = 0; i < ProcPartSize; i++) {
            for (j = 0; j < ProcPartSize; j++) {
                for (k = 0; k < dim; k++) {
                    temp += bufA[i * dim + k] * bufB[j * dim + k];
                }
                if (ProcRank - p >= 0)
                    ind = ProcRank - p;
                else ind = (ProcNum - p + ProcRank);
                bufC[i * dim + j + ind * ProcPartSize] = temp;
                temp = 0.0;
            }
        }
    }

    MPI_Gather(bufC, ProcPartElem, mpi_type, C, ProcPartElem, mpi_type, 0, MPI_COMM_WORLD);

    delete[]bufA;
    delete[]bufB;
    delete[]bufC;
}





//a solution from the internet for checking mpi
template <typename type>
void row_mul_other(type* left, type* right, type* result, size_t sz) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    MPI_Datatype mpi_type;
    if (typeid(type) == typeid(double)) mpi_type = MPI_DOUBLE;
    else if (typeid(type) == typeid(float)) mpi_type = MPI_FLOAT;
    else if (typeid(type) == typeid(int)) mpi_type = MPI_INT;
    else throw std::runtime_error("Unsupported matrix type");


    int* counts = nullptr;
    int* displs = nullptr;
    int* counts_elems = nullptr;
    int* displs_elems = nullptr;

    int local_rows;
    int total_rows = static_cast<int>(sz);


    if (rank == 0) {
        counts = new int[size];
        displs = new int[size];
        counts_elems = new int[size];
        displs_elems = new int[size];

        int offset = 0;
        for (int i = 0; i < size; ++i) {
            counts[i] = total_rows / size;
            if (i < total_rows % size) counts[i]++;

            displs[i] = offset;
            offset += counts[i];

            counts_elems[i] = counts[i] * sz;
            displs_elems[i] = displs[i] * sz;
        }
        local_rows = counts[0];
    }


    MPI_Scatter(counts, 1, MPI_INT, &local_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);


    type* local_A = new type[local_rows * sz]();
    type* local_C = new type[local_rows * sz]();
    type* B = new type[sz * sz];


    if (rank == 0) {
        std::copy(right, right + sz * sz, B);
    }
    MPI_Bcast(B, sz * sz, mpi_type, 0, MPI_COMM_WORLD);


    MPI_Scatterv(
        left,
        counts_elems,
        displs_elems,
        mpi_type,
        local_A,
        local_rows * sz,
        mpi_type,
        0,
        MPI_COMM_WORLD
    );


    for (int i = 0; i < local_rows; ++i) {
        for (int j = 0; j < sz; ++j) {
            type sum = 0;
            for (int k = 0; k < sz; ++k) {
                sum += local_A[i * sz + k] * B[k * sz + j];
            }
            local_C[i * sz + j] = sum;
        }
    }


    MPI_Gatherv(
        local_C,
        local_rows * sz,
        mpi_type,
        result,
        counts_elems,
        displs_elems,
        mpi_type,
        0,
        MPI_COMM_WORLD
    );


    delete[] local_A;
    delete[] local_C;
    delete[] B;

    if (rank == 0) {
        delete[] counts;
        delete[] displs;
        delete[] counts_elems;
        delete[] displs_elems;
    }
}