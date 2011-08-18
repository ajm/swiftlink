using namespace std;

#include <cstdio>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>
#include <string>

#include "trait.h"
#include "peel_matrix.h"


unsigned int PeelMatrix::generate_index(PeelMatrixKey& pmk) {
    unsigned int index = 0;

    for(unsigned int i = 0; i < num_keys; ++i) {
        index += (offsets[i] * pmk.get(keys[i]));
    }
    
    return index;
}

PeelMatrix::PeelMatrix(unsigned int num_dim, unsigned int val_dim) :
    num_keys(0),
    keys(NULL),
    offsets(NULL),
    number_of_dimensions(num_dim),
    values_per_dimension(val_dim),
    size((unsigned int) pow(static_cast<double>(values_per_dimension), static_cast<double>(number_of_dimensions))),
    data(NULL) {
        
//    size = (unsigned int) pow(static_cast<double>(values_per_dimension), 
//	       			static_cast<double>(number_of_dimensions));

    data = new double[size];
    
    for(unsigned i = 0; i < size; ++i) {
        data[i] = 0.0;
    }
}

PeelMatrix::PeelMatrix(const PeelMatrix& rhs) :
    num_keys(rhs.num_keys),
    keys(NULL),
    offsets(NULL),
    number_of_dimensions(rhs.number_of_dimensions),
    values_per_dimension(rhs.values_per_dimension),
    size(rhs.size),
    data(NULL) {
    
    data = new double[size];
    copy(rhs.data, rhs.data + size, data);
    
    keys = new unsigned int[rhs.num_keys];
    copy(rhs.keys, rhs.keys + num_keys, keys);
    
    offsets = new unsigned int[rhs.num_keys];
    copy(rhs.offsets, rhs.offsets + num_keys, offsets);
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
        
            delete[] offsets;
            offsets = new unsigned int[rhs.num_keys];
        }
                
        num_keys = rhs.num_keys;
        copy(rhs.keys,    rhs.keys    + num_keys, keys);
        copy(rhs.offsets, rhs.offsets + num_keys, offsets);
    }

    return *this;
}

PeelMatrix::~PeelMatrix() {
    delete[] data;
    delete[] keys;
    delete[] offsets;
}

void PeelMatrix::set_keys(vector<unsigned int>& k) {
    //keys = k;
    //sort(keys.begin(), keys.end()); // needed to do a comparison later...
    num_keys = k.size();
    keys = new unsigned int[num_keys];
    copy(k.begin(), k.end(), keys);
    sort(keys, keys + num_keys);
    init_offsets();
}

void PeelMatrix::init_offsets() {
    offsets = new unsigned int[num_keys];
    
    for(unsigned int i = 0; i < num_keys; ++i) {
        offsets[i] = (unsigned int) pow(static_cast<double>(values_per_dimension), static_cast<double>(i));
    }
}

/*
bool PeelMatrix::is_legal(PeelMatrixKey& pmk) {
    return pmk.check_keys(keys);
}
*/
double PeelMatrix::get(PeelMatrixKey& pmk) {
    return data[generate_index(pmk)];
}

void PeelMatrix::set(PeelMatrixKey& pmk, double value) {
    data[generate_index(pmk)] = value;
}

void PeelMatrix::add(PeelMatrixKey& pmk, double value) {
    data[generate_index(pmk)] += value;
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
    }
    
    for(unsigned i = 0; i < size; ++i) {
        data[i] /= matrix_sum;
    }
}

/*
// XXX stolen from Rfunction
void PeelMatrix::generate_key(PeelMatrixKey& pmatrix_index, vector<unsigned int>& assignments) {
    pmatrix_index.reassign(keys, assignments);
}

void PeelMatrix::print_keys() {
    for(unsigned i = 0; i < keys.size(); ++i)
        printf("%d ", int(keys[i]));
    printf("\n");
}
    
void PeelMatrix::print() {
    PeelMatrixKey k;
    vector<unsigned int> q;
    unsigned int ndim = keys.size();
    unsigned int tmp;
    unsigned int i;
    
    // initialise to the first element of matrix
    for(i = 0; i < ndim; ++i) {
        q.push_back(0);
    }
    
    // enumerate all elements in ndim-dimenstional matrix
    while(not q.empty()) {
            
        if(q.size() == ndim) {
            generate_key(k, q);
            
            k.print();
            printf(" (offset = %d)\t", generate_index(k));
            printf(" := %e\n", get(k));
            //printf(" := %.4f\n", get(k));
        }
            
        tmp = q.back() + 1;
        q.pop_back();
        
        if(tmp < 4) {
            q.push_back(tmp);
            tmp = ndim - q.size();
            // fill out rest with zeroes
            for(i = 0; i < tmp; ++i) {
                q.push_back(0);
            }
        }
    }
}
*/

