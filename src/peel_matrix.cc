#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>
#include <algorithm>
#include <string>

#include "trait.h"
#include "peel_matrix.h"

using namespace std;


PeelMatrix::PeelMatrix(unsigned int num_dim, unsigned int val_dim) :
    keys(),
    number_of_dimensions(num_dim),
    values_per_dimension(val_dim),
    size((unsigned int) pow(static_cast<double>(values_per_dimension), static_cast<double>(number_of_dimensions))),
    data(NULL) {
    
    data = new double[size];
    reset();
}

void PeelMatrix::reset() {
    fill(data, data + size, 0.0);
}

PeelMatrix::PeelMatrix(const PeelMatrix& rhs) :
    keys(rhs.keys),
    number_of_dimensions(rhs.number_of_dimensions),
    values_per_dimension(rhs.values_per_dimension),
    size(rhs.size),
    data(NULL) {
    
    data = new double[size];
    copy(rhs.data, rhs.data + size, data);
}

PeelMatrix& PeelMatrix::operator=(const PeelMatrix& rhs) {

    if(this != &rhs) {
        number_of_dimensions = rhs.number_of_dimensions;
        values_per_dimension = rhs.values_per_dimension;

        if(size != rhs.size) {
            delete[] data;
            data = new double[rhs.size];
        }

        size = rhs.size;
        copy(rhs.data, rhs.data + size, data);
        
        keys = rhs.keys;
    }

    return *this;
}

PeelMatrix::~PeelMatrix() {
    delete[] data;
}

void PeelMatrix::set_keys(vector<unsigned int>& k) {
    keys = k;
}

double PeelMatrix::get_result() {
    if(size != 1) {
        fprintf(stderr, "Cannot get result from an intermediate r-function\n");
        abort();
    }

    return data[0];
}

double PeelMatrix::sum() {
    double tmp = 0.0;
    
    for(unsigned i = 0; i < size; ++i) {
        tmp += data[i];
    }
        
    return tmp;
}

void PeelMatrix::normalise() {
    double matrix_sum = sum();
    
    if(matrix_sum == 0.0) {
        fprintf(stderr, "error: %s:%d, zero sum in normalise\n", __func__, __LINE__);
        abort();
    }
    
    for(unsigned i = 0; i < size; ++i) {
        data[i] /= matrix_sum;
    }
}

