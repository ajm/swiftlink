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
class GeneticMap;

class Rfunction {

    PeelMatrix pmatrix;    
    PeelOperation peel;
    unsigned int num_alleles; // could be 3 for L-sampler or 4 for peeling
    GeneticMap* map;
    Pedigree* ped;
    Person* pivot;
    
    vector<unsigned int> missing;
    vector<unsigned int> additional;
    

    void generate_key(PeelMatrixKey& pmatrix_index, vector<unsigned int>& assignments);
    bool affected_trait(enum phased_trait pt, int allele);
    enum phased_trait get_phased_trait(
                    enum phased_trait m, 
                    enum phased_trait p, 
                    int maternal_allele, 
                    int paternal_allele
                );
    double get_disease_probability(enum phased_trait pt);
    double get_recombination_probability(
                    SimwalkDescentGraph* dg, 
                    unsigned int locus_index, 
                    int maternal_allele, 
                    int paternal_allele
                );
    void evaluate_child_peel(
                    PeelMatrixKey& pmatrix_index, 
                    PeelMatrix* prev_matrix, 
                    SimwalkDescentGraph* dg, 
                    unsigned int locus_index
                );
    void evaluate_partner_peel(
                    PeelMatrixKey& pmatrix_index, 
                    PeelMatrix* prev_matrix
                );
/*
    void evaluate_last_peel(
                    PeelMatrixKey& pmatrix_index, 
                    PeelMatrix* prev_matrix
                );
*/
    void evaluate_element(
                    PeelMatrixKey& pmatrix_index, 
                    PeelMatrix* prev_matrix, 
                    SimwalkDescentGraph* dg, 
                    unsigned int locus_index
                );

 public :
    Rfunction(PeelOperation po, Pedigree* p, GeneticMap* m, unsigned int alleles);
    
    PeelMatrix* get_matrix() { return &pmatrix; }
    
    bool evaluate(PeelMatrix* previous_matrix, SimwalkDescentGraph* dg, unsigned int locus_index);
};

#endif

