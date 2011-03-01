#ifndef LKG_RFUNCTION_H_
#define LKG_RFUNCTION_H_

using namespace std;

#include <vector>

#include "peel_matrix.h"
#include "peeling.h"
#include "trait.h"


class Pedigree;
class Person;
class SimwalkDescentGraph;

class Rfunction {

    PeelMatrix pmatrix;    
    PeelOperation peel;
    unsigned int num_alleles; // could be 3 for L-sampler or 4 for peeling
    Pedigree* ped;
    Person* pivot;

    void generate_key(PeelMatrixKey& pmatrix_index, vector<unsigned int>& assignments);
    bool affected_trait(enum phased_trait pt, int allele);
    double get_disease_probability(enum phased_trait m, enum phased_trait p, int maternal_allele, int paternal_allele);
    double get_recombination_probability(SimwalkDescentGraph* dg, int maternal_allele, int paternal_allele);
    void evaluate_child_peel(PeelMatrixKey& pmatrix_index, PeelMatrix* prev_matrix, SimwalkDescentGraph* dg);
    void evaluate_element(PeelMatrixKey& pmatrix_index, PeelMatrix* prev_matrix, SimwalkDescentGraph* dg);

 public :
    Rfunction(PeelOperation po, Pedigree* p, unsigned int alleles);
    
    PeelMatrix* get_matrix() { return &pmatrix; }
    
    bool evaluate(PeelMatrix* previous_matrix, SimwalkDescentGraph* dg);
};

#endif

