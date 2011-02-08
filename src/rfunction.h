#ifndef LKG_RFUNCTION_H_
#define LKG_RFUNCTION_H_

using namespace std;

#include "peel_matrix.h"

class PeelOperation;
class Pedigree;
class Person;

class Rfunction {

    PeelMatrix pmatrix;    
    PeelOperation& peel;
    unsigned int num_alleles; // could be 3 for L-sampler or 4 for peeling
    
    Rfunction* prev; // can be null
    Pedigree* ped;
    Person* pivot;
    
    
 public :
    Rfunction(PeelOperation& p, unsigned int alleles, Rfunction* rfun, Pedigree* p) 
        : pmatrix(p.get_cutset(), alleles), peel(p), num_alleles(alleles), prev(rfun), ped(p) {
        
        pivot = ped->get_by_index(peel.get_pivot());
    }
    ~Rfunction() {}
    
    void evaluate();
};

#endif

