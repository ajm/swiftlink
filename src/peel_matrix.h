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
    
    /*
    unsigned int generate_index(vector<int>& index) const {
        unsigned int tmp = 0;
    
        for(unsigned int i = 0; i < num_keys; ++i) {
        
            if(index[keys[i]] == -1) {
                fprintf(stderr, "%d was -1\n", keys[i]);
                abort();
            }
        
            tmp += (index[keys[i]] * (1 << (2 * i)));
        }
        
        return tmp;
    }
    */
    
    //unsigned int generate_index(vector<int>& index) const;
    
    inline int generate_index(vector<int>& index) const {
        int tmp = 0;
        
        for(int i = 0; i < int(num_keys); ++i) {
            tmp += (index[keys[i]] * (1 << (2 * i)));
            //tmp += (index[keys[i]] * offsets[i]);
        }
        
        return tmp;
    }
    
    double get(vector<int>& pmk) const {
        /*
        fprintf(stderr, "PeelMatrix::get ");
        for(unsigned int i = 0; i < num_keys; ++i) {
            fprintf(stderr, "[%d %d] ", keys[i], pmk[keys[i]]);
        }
        
        int index = generate_index(pmk);
        double x = data[index];
        
        fprintf(stderr, "= %f (index=%d)\n", x, index);
        
        return x;
        */
        return data[generate_index(pmk)];
    }
    
    double get(unsigned int pmk) const {
        return data[pmk];
    }
    
    void set(unsigned int pmk, double value) {
        data[pmk] = value;
    }
    
    void add(unsigned int pmk, double value) {
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

