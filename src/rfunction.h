#ifndef LKG_RFUNCTION_H_
#define LKG_RFUNCTION_H_

using namespace std;

#include <vector>

#include "types.h"
#include "peel_matrix.h"
#include "peeling.h"
#include "genetic_map.h"


#define NUM_ALLELES 4

class Pedigree;
class Person;
class DescentGraph;


class Rfunction {

 protected :
    GeneticMap* map;
    Pedigree* ped;
    double offset;
    PeelMatrix pmatrix;
    PeelMatrix pmatrix_presum;
    PeelOperation* peel;
    Rfunction* previous_rfunction1;
    Rfunction* previous_rfunction2;
    unsigned int locus;
    double theta;
    double antitheta;
    double theta2;
    double antitheta2;
    vector<vector<int> > indices;
    unsigned int index_offset;
    unsigned int size;
    unsigned int peel_id;
    
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
    bool legal_genotype(unsigned personid, enum phased_trait g);   
    virtual double get_trait_probability(unsigned person_id, enum phased_trait pt)=0;
    virtual void evaluate_child_peel(unsigned int pmatrix_index, DescentGraph* dg)=0;
    virtual void evaluate_parent_peel(unsigned int pmatrix_index, DescentGraph* dg)=0;
    void evaluate_partner_peel(unsigned int pmatrix_index);
    void evaluate_element(unsigned int pmatrix_index, DescentGraph* dg);

 public :
    Rfunction(Pedigree* p, GeneticMap* m, unsigned int locus, PeelOperation* po, Rfunction* prev1, Rfunction* prev2);
    Rfunction(const Rfunction& r);
    Rfunction& operator=(const Rfunction& rhs);
    virtual ~Rfunction() {}
    
    double get(vector<int>& index) {
        return pmatrix.get(index);
    }
    
    void evaluate(DescentGraph* dg, double offset);
    
    void set_locus(unsigned int locus) {
        this->locus = locus;
    
        if(locus != (map->num_markers() - 1)) {
            theta = map->get_theta(locus);
            antitheta = map->get_inversetheta(locus);
        }
    
        if(locus != 0) {
            theta2 = map->get_theta(locus-1);
            antitheta2 = map->get_inversetheta(locus-1);
        }
    }
    
    double get_result() { 
        return pmatrix.get_result();
    }
};

#endif

