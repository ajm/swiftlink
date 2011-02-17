#ifndef LKG_RFUNCTION_H_
#define LKG_RFUNCTION_H_

using namespace std;

#include <deque>

#include "peel_matrix.h"
#include "peeling.h"
#include "pedigree.h"


class Person;

class Rfunction {

    PeelMatrix pmatrix;    
    unsigned int num_alleles; // could be 3 for L-sampler or 4 for peeling
    PeelOperation peel;
    Pedigree* ped;
    Person* pivot;

    void generate_key(PeelMatrixKey& pmatrix_index, deque<unsigned int>& assignments);
    void evaluate_element(PeelMatrixKey& pmatrix_index, PeelMatrix* prev_matrix);

 public :
    Rfunction(PeelOperation po, Pedigree* p, unsigned int alleles)
        : pmatrix(po.get_cutset_size(), alleles), 
        peel(po), 
        num_alleles(alleles), 
        ped(p) {
        
        pivot = ped->get_by_index(peel.get_pivot());
        pmatrix.set_keys(peel.get_cutset());
    }
    
    PeelMatrix* get_matrix() { return &pmatrix; }
    
    bool evaluate(PeelMatrix* previous_matrix);
};

#endif

