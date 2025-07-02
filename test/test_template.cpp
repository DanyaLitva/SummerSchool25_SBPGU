//UNN
//Litvyakov D. D.
//

#include "template.h"
#include <gtest.h>

const size_t MSIZE = 1024;

template <typename type>
class TestTMatrix : public ::testing::Test
{
protected:
    void SetUp()
    {
        lmatrix = new type[msize * msize];
        rmatrix = new type[msize * msize];
        resmatrix1 = new type[msize * msize];
        resmatrix2 = new type[msize * msize];
        generatematrix(lmatrix, msize);
        generatematrix(rmatrix, msize);
        basedmul(lmatrix, rmatrix, resmatrix1, msize);
    }
    void TearDown()
    {
        delete lmatrix;
        delete rmatrix;
        delete resmatrix1;
        delete resmatrix2;
    }
    type* lmatrix = nullptr;
    type* rmatrix = nullptr;
    type* resmatrix1 = nullptr;
    type* resmatrix2 = nullptr;
    size_t msize = MSIZE;
    using param = type;
};
TYPED_TEST_CASE_P(TestTMatrix);



TYPED_TEST_P(TestTMatrix, cache_friendly)
{
    int myid, numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if (myid == 0) {
        mul(this->lmatrix, this->rmatrix, this->resmatrix2, this->msize);
        EXPECT_TRUE(checkerrormatrix(this->resmatrix1, this->resmatrix2, this->msize) < 1.e-6);
    }
}

TYPED_TEST_P(TestTMatrix, decomp_mult) {
    int myid, numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    decMatrixMultiplicationMPI(this->lmatrix, this->rmatrix, this->resmatrix2, this->msize);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myid == 0) {
        EXPECT_TRUE(checkerrormatrix(this->resmatrix1, this->resmatrix2, this->msize) < 1.e-6);
    }
}

TYPED_TEST_P(TestTMatrix, strip_mult) {
    int myid, numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    stripMatrixMultiplicationMPI(this->lmatrix, this->rmatrix, this->resmatrix2, this->msize);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myid == 0) {
        EXPECT_TRUE(checkerrormatrix(this->resmatrix1, this->resmatrix2, this->msize) < 1.e-6);
    }
}






REGISTER_TYPED_TEST_CASE_P(TestTMatrix, cache_friendly, decomp_mult, strip_mult);

typedef ::testing::Types<int, float, double> MatrixTypes;
INSTANTIATE_TYPED_TEST_CASE_P(MatrixTypesInstantiation, TestTMatrix, MatrixTypes);