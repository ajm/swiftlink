#ifndef LKG_RFUNCTION_H_
#define LKG_RFUNCTION_H_

using namespace std;

#include <vector>

#include "types.h"
#include "peel_matrix.h"
#include "peeling.h"


#define NUM_ALLELES 4

class Pedigree;
class Person;
class DescentGraph;
class GeneticMap;


class Rfunction {

 protected :
    GeneticMap* map;
    Pedigree* ped;
    double offset;
    PeelMatrix pmatrix;
    PeelMatrix pmatrix_presum;
    PeelOperation peel;
    Rfunction* previous_rfunction1;
    Rfunction* previous_rfunction2;
    bool function_used;
    double theta;
    double antitheta;
    double theta2;
    double antitheta2;
        
    //enum trait get_trait(enum phased_trait p, enum parentage parent);
    inline enum trait get_trait(enum phased_trait p, enum parentage parent) {
        switch(parent) {
            case MATERNAL:
                return (((p == TRAIT_UU) or (p == TRAIT_UA)) ? TRAIT_U : TRAIT_A);    
            case PATERNAL:
                return (((p == TRAIT_UU) or (p == TRAIT_AU)) ? TRAIT_U : TRAIT_A);
        
            default:
                break;
        }
        abort();
    }
    
    void get_traits(enum phased_trait p, enum trait& mat, enum trait& pat);
    void summation(PeelMatrixKey& pmatrix_index, unsigned personid);
        
 private :   
    void generate_key(PeelMatrixKey& pmatrix_index, vector<unsigned int>& assignments);
    virtual double get_trait_probability(unsigned person_id, enum phased_trait pt, unsigned locus)=0;
    virtual void evaluate_child_peel(
                    PeelMatrixKey& pmatrix_index, 
                    DescentGraph* dg, 
                    unsigned locus)=0;
    virtual void evaluate_parent_peel(
                    PeelMatrixKey& pmatrix_index, 
                    DescentGraph* dg,
                    unsigned locus)=0;
    void evaluate_partner_peel(
                    PeelMatrixKey& pmatrix_index, 
                    unsigned locus);
    void evaluate_element(
                    PeelMatrixKey& pmatrix_index, 
                    DescentGraph* dg, 
                    unsigned int locus_index);

 public :
    Rfunction(PeelOperation po, Pedigree* p, GeneticMap* m, Rfunction* prev1, Rfunction* prev2);
    Rfunction(const Rfunction& r);
    Rfunction& operator=(const Rfunction& rhs);
    virtual ~Rfunction() {}
    
    //PeelMatrix* get_matrix() { return &pmatrix; }
    double get(PeelMatrixKey& pmk) { return pmatrix.get(pmk); }
    
    void evaluate(DescentGraph* dg, unsigned locus, double offset);

    //void print() { pmatrix_presum.print(); }
    //void print_keys() { pmatrix.print_keys(); }
    
    bool is_used();
    void set_used();
    bool contains_node(unsigned node);
    bool contains_cutnodes(vector<unsigned>& nodes);    
};

#endif

