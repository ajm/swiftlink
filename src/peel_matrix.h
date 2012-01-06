#ifndef LKG_PEELMATRIX_H_
#define LKG_PEELMATRIX_H_

using namespace std;

#include <cstdio>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>
#include <string>

#include "trait.h"

/*
class PeelMatrixKey {
    
    unsigned int num_keys;
    enum phased_trait* key;
    
 public :
    PeelMatrixKey(unsigned max_keys) :
        num_keys(max_keys), 
        key(NULL) {
    
        key = new enum phased_trait[num_keys];
        for(unsigned i = 0; i < num_keys; ++i)
            key[i] = TRAIT_UU;
    }
    
    ~PeelMatrixKey() {
        delete[] key;
    }

    PeelMatrixKey(const PeelMatrixKey& rhs) :
        num_keys(rhs.num_keys), 
        key(NULL) {
        
        key = new enum phased_trait[num_keys];
        copy(rhs.key, rhs.key + num_keys, key);
    }

    PeelMatrixKey& operator=(const PeelMatrixKey& rhs) {
        if(this != &rhs) {
            if(rhs.num_keys != num_keys) {
                num_keys = rhs.num_keys;
                delete[] key;
                key = new enum phased_trait[num_keys];
            }
            
            copy(rhs.key, rhs.key + num_keys, key);
        }

        return *this;
    }
    
    void reassign(vector<unsigned int>& cutset, vector<unsigned int>& assignments) {
        for(unsigned int i = 0; i < cutset.size(); ++i) {
            add(
                cutset[i], 
                static_cast<enum phased_trait>(assignments[i])
            );
        }
    }

    void add(unsigned int k, enum phased_trait value) {
        key[k] = value;
    }

    inline enum phased_trait get(unsigned int i) {
        return key[i];
    }

    void print() {
        
    }
    
    void raw_print() {
        for(unsigned int i = 0; i < num_keys; ++i)
            printf("%d ", (int) key[i]);
        printf("\n");
    }
};
*/
class PeelMatrix {

    unsigned int num_keys;
    unsigned int* keys;
    //unsigned int* offsets;
    unsigned int number_of_dimensions;
    unsigned int values_per_dimension;
    unsigned int size;
    double* data;
    
    void init_offsets();
    
 public :
    PeelMatrix(unsigned int num_dim, unsigned int val_dim);
    PeelMatrix(const PeelMatrix& rhs);
    PeelMatrix& operator=(const PeelMatrix& rhs);
    ~PeelMatrix();
    
    
    void set_keys(vector<unsigned int>& k);

    double get_result();
    double sum();
    void normalise();
    
    unsigned int generate_index(vector<int>& index) const {
        unsigned int tmp = 0;
    
        for(unsigned int i = 0; i < num_keys; ++i) {
            tmp += (index[keys[i]] * (1 << (2 * i)));
        }
        
        return tmp;
    }
    
    double get(vector<int>& pmk) const {
        return data[generate_index(pmk)];
    }
    
    double get(unsigned int pmk) const {
        return data[pmk];
    }
    
    void set(unsigned int pmk, double value) {
        //if(value != 0.0)
            data[pmk] = value;
    }
    
    void add(unsigned int pmk, double value) {
        //if(value != 0.0)
            data[pmk] += value;
    }
    
    void reset();
    
    void raw_print() {
        for(unsigned int i = 0; i < size; ++i) {
            printf("%.3f\n", data[i]);
        }
        printf("\n");
    }
};

#endif

