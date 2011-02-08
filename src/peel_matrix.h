#ifndef LKG_PEELMATRIX_H_
#define LKG_PEELMATRIX_H_

using namespace std;

#include <cmath>
#include <map>

#include "genotypes.h"


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
        if(*this != rhs) {
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
    vector<unsigned int>& keys;
    vector<unsigned int> offsets;
    unsigned int number_of_dimensions;
    unsigned int values_per_dimension;
    unsigned int size;
    double* data;
    
    unsigned int generate_index(PeelMatrixKey& pmk) {
        unsigned int index = 0;

        for(unsigned int i = 0; i < keys.size(); ++i) {
            index += (offset[i] * pmk.get(keys[i]));
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
    PeelMatrix(vector<unsigned int>& all_keys, unsigned int val_dim)
        : keys(all_keys), 
        number_of_dimensions(all_keys.size()), 
        values_per_dimension(val_dim) {
        
        init_offsets();
        size = pow(values_per_dimension, number_of_dimensions);
        data = new double[size];
    }
    
    ~PeelMatrix() {
        delete[] data;
    }

    bool is_legal(PeelMatrixKey& pmk) {
        return pmk.check_keys();
    }

    double get(PeelMatrixKey& pmk) {
        return data[generate_index(pmk)];
    }

    void set(PeelMatrixKey& pmk, double value) {
        data[generate_index(pmk)] = value;
    }
    
    // TODO overload '[]' operator?
};

#endif

