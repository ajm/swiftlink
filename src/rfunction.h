#ifndef LKG_RFUNCTION_H_
#define LKG_RFUNCTION_H_

using namespace std;

#include "peel_matrix.h"


class PeelOperation;
class Pedigree;
class Person;

class Rfunction {

    PeelMatrix pmatrix;    
    unsigned int num_alleles; // could be 3 for L-sampler or 4 for peeling
    PeelOperation peel;
    Pedigree* ped;
    Person* pivot;
    
 public :
    Rfunction(PeelOperation p, Pedigree* p, unsigned int alleles)
        : pmatrix(p.get_cutset(), alleles), peel(p), num_alleles(alleles), ped(p) {
        
        pivot = ped->get_by_index(peel.get_pivot());
    }
    
    PeelMatrix* get_matrix() { return &pmatrix; }
    
    bool evaluate(PeelMatrix* previous_matrix);
};

#endif

