//UNN
//Litvyakov D. D.
//

#pragma once
#include <random>

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
        for (size_t k = 0; k < sz; ++k) {  
            type tempval = lefttemp[k];
            type* resulttemp = &result[i * sz];
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