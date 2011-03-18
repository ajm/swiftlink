#ifndef LKG_PEELMATRIX_H_
#define LKG_PEELMATRIX_H_

using namespace std;

#include <cstdio>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>
#include <string>

//#include "genotype.h"
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
            offsets.push_back(pow(static_cast<double>(values_per_dimension), 
                                  static_cast<int>(i)));
        }
    }
    
 public :
    PeelMatrix(unsigned int num_dim, unsigned int val_dim)
        : /*keys(all_keys), */
        /*number_of_dimensions(all_keys.size()), */
        number_of_dimensions(num_dim),
        values_per_dimension(val_dim) {
        
        //init_offsets();
        size = pow(static_cast<double>(values_per_dimension), 
                   static_cast<int>(number_of_dimensions));
        data = new double[size];
        
        //printf("size = %d\n", size);
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
        delete[] data;
    }

    bool key_intersection(PeelMatrix* pm, 
        vector<unsigned int>& missing, vector<unsigned int>& additional) {
        
        unsigned int ikey;

        missing = keys;

        if(!pm) {
            return false;
        }

        additional = pm->keys;
/*
        for(unsigned int i = 0; i < missing.size(); ++i) {
            fprintf(stderr, "missing[%d] = %d\n", i, missing[i]);
        }
        for(unsigned int i = 0; i < additional.size(); ++i) {
            fprintf(stderr, "additional[%d] = %d\n", i, additional[i]);
        }
*/
        // this looks bad (n^2), but the number of dimensions is pretty 
        // constrained
        for(unsigned int i = 0; i < keys.size(); ++i) {
            ikey = keys[i];

            if(binary_search(additional.begin(), additional.end(), ikey)) {
                missing.erase(find(missing.begin(), missing.end(), ikey));
                additional.erase(find(additional.begin(), additional.end(), ikey));
            }
        }
/*
        for(unsigned int i = 0; i < missing.size(); ++i) {
            fprintf(stderr, "missing[%d] = %d\n", i, missing[i]);
        }
        for(unsigned int i = 0; i < additional.size(); ++i) {
            fprintf(stderr, "additional[%d] = %d\n", i, additional[i]);
        }
*/
        return (missing.size() == 0) and (additional.size() == 0);
    }

    void set_keys(vector<unsigned int>& k) {
        keys = k;
        sort(keys.begin(), keys.end()); // needed to do a comparison later...
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
/*
    void set_raw(unsigned int index, double value) {
        data[index] = value;
    }
    
    double get_raw(unsigned int index) {
        return data[index];
    }
*/
    double sum() {
        double tmp = 0.0;
        
        for(unsigned i = 0; i < size; ++i) {
            tmp += data[i];
        }
        
        return tmp;
    }

    // XXX stolen from Rfunction
    void generate_key(PeelMatrixKey& pmatrix_index, vector<unsigned int>& assignments) {
        pmatrix_index.reassign(keys, assignments);
    }
    
    void print() {
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
                printf(" := %f\n", get(k));
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
};

#endif

