#ifndef LKG_RFUNCTION_H_
#define LKG_RFUNCTION_H_

using namespace std;

#include <cmath>
#include <deque>


class Rfunction {
    
    PeelOperation& peel;
    unsigned int num_alleles; // could be 3 for L-sampler or 4 for peeling
    double* rfunc;
    Rfunction* prev; // can be null
    Pedigree* ped;
    Person* pivot;
    
    void normalise(); // specific dimension(s)?
    
    
 public :
    Rfunction(PeelOperation& p, unsigned int alleles, Rfunction* rfun, Pedigree* p) 
        : peel(p), num_alleles(alleles), prev(rfun), ped(p) {
        rfunc = new double[static_cast<int>(pow(num_alleles, peel.get_cutset_size()))];
        pivot = ped->get_by_index(peel.get_pivot());
    }
    
    ~Rfunction() {
        delete[] rfunc;
    }
    
    void evaluate();
};

#endif

