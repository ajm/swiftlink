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

 protected :
    GeneticMap* map;
    Pedigree* ped;
    double offset;
    
 private :
    PeelMatrix pmatrix;
    PeelMatrix pmatrix_presum;
    
    PeelOperation peel;
    
    
    Rfunction* previous_rfunction1;
    Rfunction* previous_rfunction2;
    bool function_used;
    unsigned function_index;
    
    
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
    virtual double get_trait_probability(unsigned person_id, enum phased_trait pt, unsigned locus)=0;
    virtual double get_recombination_probability(
                    DescentGraph* dg, 
                    unsigned locus_index,
                    unsigned person_id,
                    int maternal_allele, 
                    int paternal_allele
                )=0;
    double get_recombination_probability_between_markers(
                    DescentGraph* dg, 
                    unsigned locus_index,
                    unsigned person_id,
                    int maternal_allele, 
                    int paternal_allele
                );
    void summation(PeelMatrixKey& pmatrix_index, unsigned personid);
    void evaluate_child_peel(
                    PeelMatrixKey& pmatrix_index, 
                    DescentGraph* dg, 
                    unsigned locus);
    void evaluate_parent_peel(
                    PeelMatrixKey& pmatrix_index, 
                    DescentGraph* dg,
                    unsigned locus);
    void evaluate_partner_peel(
                    PeelMatrixKey& pmatrix_index, 
                    unsigned locus);
    void evaluate_element(
                    PeelMatrixKey& pmatrix_index, 
                    DescentGraph* dg, 
                    unsigned int locus_index);

 public :
    Rfunction(PeelOperation po, Pedigree* p, GeneticMap* m,
        vector<Rfunction*>& previous_functions, unsigned index);
                
    virtual ~Rfunction() {}
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

