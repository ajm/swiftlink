#ifndef LKG_RFUNCTION_H_
#define LKG_RFUNCTION_H_

using namespace std;

#include <cmath>


class Rfunction {

    PeelOperation& peel;
    unsigned int num_alleles;
    double* rfunc;

    void normalise() {} // specific dimension(s)?
    
    
 public :
    Rfunction(PeelOperation& p, unsigned int alleles) 
        : peel(p), num_alleles(alleles) {
        rfunc = new double[static_cast<int>(pow(num_alleles, peel.get_cutset_size()))];
    }

    ~Rfunction() {
        delete[] rfunc;
    }

    void evaluate() {
        
    }
};

#endif

