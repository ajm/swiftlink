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


//template<class A, class B>
class PeelMatrixKey {

    // the key to the peel matrix itself is a map that needs to be 
    // transformed into an unsigned integer
    // nodeid,unphased_genotype
    // nodeid,phased_genotype
//    map<A,B> key;
    map<unsigned int, enum phased_trait> key;
    
 public :
    PeelMatrixKey() {}
    
    PeelMatrixKey(vector<unsigned int>& cutset, vector<unsigned int>& assignments) {
        reassign(cutset, assignments);
    }

    ~PeelMatrixKey() {}

    PeelMatrixKey(const PeelMatrixKey& pmk) {
        key = pmk.key;
    }

    PeelMatrixKey& operator=(const PeelMatrixKey& rhs) {
        if(this != &rhs) {
            key = rhs.key;
        }

        return *this;
    }

//    void add(A& key, B& val);
//    bool remove(A& key);

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
    
    void remove(unsigned int k) {
        key.erase(k);
    }

    enum phased_trait get(unsigned int i) {
        return key[i];
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

    void print() {
        map<unsigned int, enum phased_trait>::iterator it;
        
        for(it = key.begin(); it != key.end(); it++) {
            printf("%d=%d ", (*it).first, (*it).second);
        }
    }
};

class PeelMatrix {

    vector<unsigned int> keys;
    vector<unsigned int> offsets;
    unsigned int number_of_dimensions;
    unsigned int values_per_dimension;
    unsigned int size;
    double* data;
    
    unsigned int generate_index(PeelMatrixKey& pmk);
    void init_offsets();
    
 public :
    PeelMatrix(unsigned int num_dim, unsigned int val_dim);
    PeelMatrix(const PeelMatrix& rhs);
    PeelMatrix& operator=(PeelMatrix& rhs);
    ~PeelMatrix();

    bool key_intersection(
            PeelMatrix* pm, 
            vector<unsigned int>& missing, 
            vector<unsigned int>& additional
        );
    void set_keys(vector<unsigned int>& k);
    bool is_legal(PeelMatrixKey& pmk);
    double get(PeelMatrixKey& pmk);
    void set(PeelMatrixKey& pmk, double value);
    double get_result();
    
    void generate_key(PeelMatrixKey& pmatrix_index, vector<unsigned int>& assignments);
    void print();};

#endif

