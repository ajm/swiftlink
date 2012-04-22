using namespace std;

#include <cstdio>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>
#include <algorithm>
#include <string>

#include "trait.h"
#include "peel_matrix.h"


PeelMatrix::PeelMatrix(unsigned int num_dim, unsigned int val_dim) :
    num_keys(0),
    keys(NULL),
    number_of_dimensions(num_dim),
    values_per_dimension(val_dim),
    size((unsigned int) pow(static_cast<double>(values_per_dimension), static_cast<double>(number_of_dimensions))),
    data(NULL) {
    
    data = new double[size];
    reset();
}

void PeelMatrix::reset() {
    /*
    for(unsigned i = 0; i < size; ++i) {
        data[i] = 0.0;
    }
    
    memset(data, 0, sizeof(double) * size);
    */
    fill(data, data + size, 0.0);
}

PeelMatrix::PeelMatrix(const PeelMatrix& rhs) :
    num_keys(rhs.num_keys),
    keys(NULL),
    number_of_dimensions(rhs.number_of_dimensions),
    values_per_dimension(rhs.values_per_dimension),
    size(rhs.size),
    data(NULL) {
    
    data = new double[size];
    copy(rhs.data, rhs.data + size, data);
    
    keys = new unsigned int[num_keys];
    copy(rhs.keys, rhs.keys + num_keys, keys);
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
        
        
        if(num_keys != rhs.num_keys) {
            delete[] keys;
            keys = new unsigned int[rhs.num_keys];
        }
                
        num_keys = rhs.num_keys;
        copy(rhs.keys, rhs.keys + num_keys, keys);
    }

    return *this;
}

PeelMatrix::~PeelMatrix() {
    delete[] data;
    delete[] keys;
}

void PeelMatrix::set_keys(vector<unsigned int>& k) {
    num_keys = k.size();
    keys = new unsigned int[num_keys];
    copy(k.begin(), k.end(), keys);
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

unsigned int PeelMatrix::generate_index(vector<int>& index) const {
    unsigned int tmp = 0;
    
    for(unsigned int i = 0; i < num_keys; ++i) {
        /*
        if(index[keys[i]] == -1) {
            fprintf(stderr, "%d was -1\n", keys[i]);
            abort();
        }
        */
        
        tmp += (index[keys[i]] * (1 << (2 * i)));
    }
    
    return tmp;
}

