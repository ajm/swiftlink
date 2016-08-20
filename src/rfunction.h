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
    unsigned int offset;
    PeelMatrix pmatrix;
    PeelMatrix pmatrix_presum;
    PeelOperation* peel;
    vector<Rfunction*> previous_rfunctions;
    unsigned int locus;
    vector<vector<int> > indices;
    vector<int>* valid_indices;
    vector<int>* valid_lod_indices;
    unsigned int index_offset;
    unsigned int size;
    unsigned int peel_id;
    double theta;
    double antitheta;
    double theta2;
    double antitheta2;
    
    double trait_cache[4];

    bool sex_linked;
    
    //enum trait get_trait(enum phased_trait p, enum parentage parent);
    //bool affected_trait(enum phased_trait pt, int allele);
    enum phased_trait get_phased_trait(enum phased_trait m, enum phased_trait p, int maternal_allele, int paternal_allele, enum sex child_sex);
    
    void normalise(double* p);
    
    inline enum trait get_trait(enum phased_trait p, enum parentage parent) {
        /*
        switch(parent) {
            case MATERNAL:
                return (((p == TRAIT_UU) or (p == TRAIT_UA)) ? TRAIT_U : TRAIT_A);    
            case PATERNAL:
                return (((p == TRAIT_UU) or (p == TRAIT_AU)) ? TRAIT_U : TRAIT_A);
            default:
                break;
        }
        
        abort();
        */
        
        switch(p) {
            case TRAIT_UU :
                return TRAIT_U;
            case TRAIT_AU :
                return (parent == MATERNAL) ? TRAIT_A : TRAIT_U;
            case TRAIT_UA :
                return (parent == MATERNAL) ? TRAIT_U : TRAIT_A;
            case TRAIT_AA :
                return TRAIT_A;
        }
    }
    
    inline bool affected_trait(enum phased_trait pt, int allele) {
        /*
        switch(allele) {
            case 0 :
                return (pt == TRAIT_AU) or (pt == TRAIT_AA);
            case 1 :
                return (pt == TRAIT_UA) or (pt == TRAIT_AA);
            default :
                break;
        }
        
        abort();
        */
        
        switch(pt) {
            case TRAIT_UU :
                return false;
            case TRAIT_AU :
                return allele == 0;
            case TRAIT_UA :
                return allele == 1;
            case TRAIT_AA :
                return true;
        }
    }
    
        
 private :
    bool legal_genotype(unsigned personid, enum phased_trait g);
    virtual void preevaluate_init(DescentGraph* dg)=0;
    virtual double get_trait_probability(unsigned person_id, enum phased_trait pt)=0;
    virtual void evaluate_child_peel(unsigned int pmatrix_index, DescentGraph* dg)=0;
    virtual void evaluate_parent_peel(unsigned int pmatrix_index, DescentGraph* dg)=0;
    void evaluate_partner_peel(unsigned int pmatrix_index);
    void evaluate_element(unsigned int pmatrix_index, DescentGraph* dg);

 public :
    Rfunction(Pedigree* p, GeneticMap* m, unsigned int locus, PeelOperation* po, vector<Rfunction*> previous, bool sex_linked);
    Rfunction(const Rfunction& r);
    Rfunction& operator=(const Rfunction& rhs);
    virtual ~Rfunction() {}
    
    double get(vector<int>& index) {
        return pmatrix.get(index);
    }
    
    void evaluate(DescentGraph* dg, unsigned int offset);
    
    double get_result() { 
        return pmatrix.get_result();
    }
};

#endif

