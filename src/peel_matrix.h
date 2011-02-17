#ifndef LKG_PEELMATRIX_H_
#define LKG_PEELMATRIX_H_

using namespace std;

#include <cstdio>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>

#include "genotype.h"


//template<class A, class B>
class PeelMatrixKey {

    // the key to the peel matrix itself is a map that needs to be 
    // transformed into an unsigned integer
    // nodeid,unphased_genotype
    // nodeid,phased_genotype
//    map<A,B> key;
    map<unsigned int, enum phased_genotype> key;
    
 public :
    PeelMatrixKey() {}
    ~PeelMatrixKey() {}

    PeelMatrixKey(const PeelMatrixKey& pmk) {
        key = pmk.key;
    }

    PeelMatrixKey& operator=(const PeelMatrixKey& rhs) {
        if(this != &rhs) {
            key = rhs.key;
        }
    }

//    void add(A& key, B& val);
//    bool remove(A& key);

    void add(unsigned int k, enum phased_genotype value) {
        key[k] = value;
    }
    
    void remove(unsigned int k) {
        key.erase(k);
    }

    int get(unsigned int i) {
        return static_cast<int>(key[i]);
    }
    
    // ensure this key can address everything for everything in the
    // vector 'keys'
    bool check_keys(vector<unsigned int>& keys) {
        
        if(keys.size() != key.size()) {
            return false;
        }

        for(unsigned int i = 0; i < keys.size(); ++i) {
            if(key.count(keys[i]) == 0) {
                return false;
            }
        }

        return true;
    }
};

class PeelMatrix {
    vector<unsigned int> keys;
    vector<unsigned int> offsets;
    unsigned int number_of_dimensions;
    unsigned int values_per_dimension;
    unsigned int size;
    double* data;
    
    unsigned int generate_index(PeelMatrixKey& pmk) {
        unsigned int index = 0;
        unsigned int tmp;

        for(unsigned int i = 0; i < keys.size(); ++i) {
            tmp = keys[i];
            //fprintf(stderr, "%d --> %d\n", i, tmp);
            index += (offsets[i] * pmk.get(keys[i]));
        }

        return index;
    }
    
    void init_offsets() {
        offsets.resize(keys.size());
        
        for(unsigned int i = 0; i < keys.size(); ++i) {
            offsets.push_back(pow(values_per_dimension, i));
        }
    }
    
 public :
    PeelMatrix(unsigned int num_dim, unsigned int val_dim)
        : /*keys(all_keys), */
        /*number_of_dimensions(all_keys.size()), */
        number_of_dimensions(num_dim),
        values_per_dimension(val_dim) {
        
        //init_offsets();
        size = pow(values_per_dimension, number_of_dimensions);
        data = new double[size];
    }

    PeelMatrix(const PeelMatrix& rhs) {
        keys = rhs.keys;
        offsets = rhs.offsets;
        number_of_dimensions = rhs.number_of_dimensions;
        values_per_dimension = rhs.values_per_dimension;
        size = rhs.size;
        data = new double[size];
        copy(rhs.data, rhs.data + size, data);
    }

    PeelMatrix& operator=(PeelMatrix& rhs) {
        
        if(this != &rhs) {
            keys = rhs.keys;
            offsets = rhs.offsets;
            number_of_dimensions = rhs.number_of_dimensions;
            values_per_dimension = rhs.values_per_dimension;

            if(size != rhs.size) {
                delete[] data;
                data = new double[rhs.size];
            }

            size = rhs.size;
            copy(rhs.data, rhs.data + size, data);
        }

        return *this;
    }

    ~PeelMatrix() {
        delete[] data; // freed twice due to shallow default copy constructor
    }

    void set_keys(vector<unsigned int>* k) {
        keys = *k;
        init_offsets();
    }

    bool is_legal(PeelMatrixKey& pmk) {
        return pmk.check_keys(keys);
    }

    double get(PeelMatrixKey& pmk) {
        return data[generate_index(pmk)];
    }

    void set(PeelMatrixKey& pmk, double value) {
        data[generate_index(pmk)] = value;
    }
    
    // TODO overload '[]' operator ?
};

#endif

