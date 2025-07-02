//UNN
//Litvyakov D. D.
//

#pragma once
#include <random>
#include "mpi.h"
#include <iostream>
#include <vector>
#include <stdexcept>


const int MinVal = -1000;
const int MaxVal = 1000;


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
type checkerrormatrix(type* left, type* right, size_t sz) {
    type temp = type(0);
    for (size_t i = 0; i < (sz * sz); ++i)
        if (std::abs(left[i] - right[i]) > temp)  temp = std::abs(left[i] - right[i]);
    return temp;
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
void decMatrixMultiplicationMPI(type*& A, type*& B, type*& C, size_t& Size) {

    int ProcRank, ProcNum;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);

    MPI_Datatype mpi_type;
    if (typeid(type) == typeid(double)) mpi_type = MPI_DOUBLE;
    else if (typeid(type) == typeid(float)) mpi_type = MPI_FLOAT;
    else if (typeid(type) == typeid(int)) mpi_type = MPI_INT;
    else throw std::runtime_error("Unsupported matrix type");

    int dim = Size;
    int rows_per_proc = dim / ProcNum;
    int remainder = dim % ProcNum;
    int my_rows = rows_per_proc;
    if (ProcRank < remainder) {
        my_rows++;
    }

    int* counts = nullptr;
    int* displs = nullptr;

    if (ProcRank == 0) {
        counts = new int[ProcNum];
        displs = new int[ProcNum];
        int offset = 0;
        for (int i = 0; i < ProcNum; i++) {
            int rows_i = rows_per_proc;
            if (i < remainder) {
                rows_i++;
            }
            counts[i] = rows_i * dim;
            displs[i] = offset;
            offset += counts[i];
        }
    }

    type* localA = new type[my_rows * dim];
    MPI_Scatterv(A, counts, displs, mpi_type, localA, my_rows * dim, mpi_type, 0, MPI_COMM_WORLD);

    type* localB = nullptr;
    if (ProcRank == 0) {
        localB = B;
    }
    else {
        localB = new type[dim * dim];
    }
    MPI_Bcast(localB, dim * dim, mpi_type, 0, MPI_COMM_WORLD);

    type* localC = new type[my_rows * dim];
    for (int i = 0; i < my_rows; i++) {
        for (int j = 0; j < dim; j++) {
            type temp = 0;
            for (int k = 0; k < dim; k++) {
                temp += localA[i * dim + k] * localB[k * dim + j];
            }
            localC[i * dim + j] = temp;
        }
    }

    MPI_Gatherv(localC, my_rows * dim, mpi_type, C, counts, displs, mpi_type, 0, MPI_COMM_WORLD);

    delete[] localA;
    delete[] localC;
    if (ProcRank != 0) {
        delete[] localB;
    }
    if (ProcRank == 0) {
        delete[] counts;
        delete[] displs;
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


template <typename type>
void stripMatrixMultiplicationMPI(type*& A, type*& B, type*& C, size_t& Size) {

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

    temp = type(0);
    for (i = 0; i < ProcPartSize; i++) {
        for (j = 0; j < ProcPartSize; j++) {
            for (k = 0; k < dim; k++)
                temp += bufA[i * dim + k] * bufB[j * dim + k];
            bufC[i * dim + j + ProcPartSize * ProcRank] = temp;
            temp = type(0);
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
        temp = type(0);
        for (i = 0; i < ProcPartSize; i++) {
            for (j = 0; j < ProcPartSize; j++) {
                for (k = 0; k < dim; k++) {
                    temp += bufA[i * dim + k] * bufB[j * dim + k];
                }
                if (ProcRank - p >= 0)
                    ind = ProcRank - p;
                else ind = (ProcNum - p + ProcRank);
                bufC[i * dim + j + ind * ProcPartSize] = temp;
                temp = type(0);
            }
        }
    }

    MPI_Gather(bufC, ProcPartElem, mpi_type, C, ProcPartElem, mpi_type, 0, MPI_COMM_WORLD);

    if (ProcRank == 0) {
        Flip(B, Size);
    }

    delete[]bufA;
    delete[]bufB;
    delete[]bufC;
}





