//UNN
//Litvyakov D. D.
//

#pragma once
#include <random>

int MinVal = -1000;
int MaxVal = 1000;


template <typename type>
void mul(type* left, type* right, type* result, size_t sz) {
    for (size_t i = 0; i < sz * sz; ++i) result[i] = type(0);

    for (size_t i = 0; i < sz; ++i) {
        for (size_t k = 0; k < sz; ++k) {  
            type temp = left[i * sz + k];
            for (size_t j = 0; j < sz; ++j) {
                result[i * sz + j] += temp * right[k * sz + j];
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