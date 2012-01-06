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
    vector<vector<int> >* indices;
    unsigned int size;
        
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
        
 private :
    bool legal_genotype(unsigned personid, unsigned locus, enum phased_trait g);   
    virtual double get_trait_probability(unsigned person_id, enum phased_trait pt, unsigned locus)=0;
    virtual void evaluate_child_peel(
                    unsigned int pmatrix_index, 
                    DescentGraph* dg, 
                    unsigned locus)=0;
    virtual void evaluate_parent_peel(
                    unsigned int pmatrix_index, 
                    DescentGraph* dg,
                    unsigned int locus)=0;
    void evaluate_partner_peel(
                    unsigned int pmatrix_index, 
                    unsigned int locus);
    void evaluate_element(
                    unsigned int pmatrix_index, 
                    DescentGraph* dg, 
                    unsigned int locus_index);

 public :
    Rfunction(PeelOperation po, Pedigree* p, GeneticMap* m, Rfunction* prev1, Rfunction* prev2);
    Rfunction(const Rfunction& r);
    Rfunction& operator=(const Rfunction& rhs);
    virtual ~Rfunction() {}
    
    double get(vector<int>& index) {
        return pmatrix.get(index);
    }
    
    void evaluate(DescentGraph* dg, unsigned locus, double offset);   
};

#endif

