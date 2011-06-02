#ifndef LKG_RFUNCTION_H_
#define LKG_RFUNCTION_H_

using namespace std;

#include <vector>

#include "peel_matrix.h"
#include "peeling.h"
#include "trait.h"


class Pedigree;
class Person;
class DescentGraph;
class GeneticMap;

class Rfunction {

    PeelMatrix pmatrix;    
    PeelOperation peel;
    unsigned num_alleles; // could be 3 for L-sampler or 4 for peeling
    GeneticMap* map;
    Pedigree* ped;
    
    // these have different meanings for different peeling operations
    //
    // CHILD_PEEL
    //      1. same as what is being peel to (or NULL, if it does not exist)
    //      2. cutset containing the child (or NULL, if child is a leaf)
    //
    // PARTNER_PEEL
    //      1. cutset containing marriage (CANNOT BE NULL)
    //      2. cutset containing the person being peeled (or NULL, if person is a founder)
    //
    // PARENT_PEEL
    //      1. cutset containing mother (or NULL, if mother is a founder)
    //      2. cutset containing father (or NULL, if father is a founder)
    //      3. cutset containing mother x father peeled upwards (or NULL if that has not happened yet)
    //
    // LAST_PEEL
    //      1. cutset containing person being peeled
    //      2. NULL
    //      
    Rfunction* previous_rfunction1;
    Rfunction* previous_rfunction2;
    Rfunction* previous_rfunction3; // only used by PARENT_PEEL
    bool function_used;
    unsigned function_index;
    
    
    bool is_used();
    void set_used();
    void find_previous_functions(vector<Rfunction*>& functions);
    bool contains_cutnodes(vector<unsigned>& nodes);
    void find_function_containing(vector<Rfunction*>& functions, 
                                  vector<unsigned>& nodes, 
                                  Rfunction** func);
    void find_child_functions(vector<Rfunction*>& functions);
    void find_partner_functions(vector<Rfunction*>& functions);
    void find_parent_functions(vector<Rfunction*>& functions);
    void find_last_functions(vector<Rfunction*>& functions);

    void generate_key(PeelMatrixKey& pmatrix_index, vector<unsigned int>& assignments);
    bool affected_trait(enum phased_trait pt, int allele);
    enum phased_trait get_phased_trait(
                    enum phased_trait m, 
                    enum phased_trait p, 
                    int maternal_allele, 
                    int paternal_allele
                );
    double get_disease_probability(unsigned person_id, enum phased_trait pt);
    double get_recombination_probability(
                    DescentGraph* dg, 
                    unsigned int locus_index,
                    unsigned int person_id,
                    int maternal_allele, 
                    int paternal_allele
                );
    
    void evaluate_child_peel(
                    PeelMatrixKey& pmatrix_index, 
                    DescentGraph* dg, 
                    unsigned int locus_index);
    void evaluate_parent_peel(
                    PeelMatrixKey& pmatrix_index, 
                    DescentGraph* dg,
                    unsigned int locus_index);
    void evaluate_partner_peel(PeelMatrixKey& pmatrix_index);
    void evaluate_last_peel(PeelMatrixKey& pmatrix_index);
    void evaluate_element(
                    PeelMatrixKey& pmatrix_index, 
                    DescentGraph* dg, 
                    unsigned int locus_index);

 public :
    Rfunction(PeelOperation po, Pedigree* p, GeneticMap* m, unsigned alleles, 
                vector<Rfunction*>& previous_functions, unsigned index);
    
    PeelMatrix* get_matrix() { return &pmatrix; }
    double get(PeelMatrixKey& pmk) { return pmatrix.get(pmk); }
    double get_result() { return pmatrix.get_result(); }
    void print() { pmatrix.print(); }
    void print_keys() { pmatrix.print_keys(); }
    void evaluate(DescentGraph* dg, unsigned int locus_index);
};

#endif

