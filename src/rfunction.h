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


enum trait_type { 
    MEIOSIS_INDICATORS, 
    DISEASE_TRAIT 
};

class Rfunction {

    PeelMatrix pmatrix;    
    PeelOperation peel;
    unsigned num_alleles; // could be 3 for L-sampler or 4 for peeling
    GeneticMap* map;
    Pedigree* ped;
    
    Rfunction* previous_rfunction1;
    Rfunction* previous_rfunction2;
    bool function_used;
    unsigned function_index;
    enum trait_type type;
    
    
    bool is_used();
    void set_used();
    void find_previous_functions(vector<Rfunction*>& functions);
    bool contains_node(unsigned node);
    bool contains_cutnodes(vector<unsigned>& nodes);
    void find_function_containing(vector<Rfunction*>& functions, 
                                  vector<unsigned>& nodes, 
                                  Rfunction** func);
    void find_child_functions(vector<Rfunction*>& functions);
    void find_generic_functions(vector<Rfunction*>& functions);

    void generate_key(PeelMatrixKey& pmatrix_index, vector<unsigned int>& assignments);
    bool affected_trait(enum phased_trait pt, int allele);
    enum phased_trait get_phased_trait(
                    enum phased_trait m, 
                    enum phased_trait p, 
                    int maternal_allele, 
                    int paternal_allele
                );
    double get_disease_probability(unsigned person_id, enum phased_trait pt);
    double get_meiosis_probability(unsigned person_id, enum phased_trait pt);
    double get_trait_probability(unsigned person_id, enum phased_trait pt);
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
                
    ~Rfunction() {}
    Rfunction(const Rfunction& r);
    Rfunction& operator=(const Rfunction& rhs);
    
    PeelMatrix* get_matrix() { return &pmatrix; }
    double get(PeelMatrixKey& pmk) { return pmatrix.get(pmk); }
    double get_result() { return pmatrix.get_result(); }
    
    void evaluate(DescentGraph* dg, unsigned locus, double offset);

    void print() { pmatrix.print(); }
    void print_keys() { pmatrix.print_keys(); }
};

#endif

