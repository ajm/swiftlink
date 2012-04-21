#ifndef LKG_TRAITRFUNCTION_H_
#define LKG_TRAITRFUNCTION_H_

#include <vector>
#include "rfunction.h"


class TraitRfunction : public Rfunction {
    
    double half_theta;
    double half_inversetheta;
    
    double get_recombination_probability(DescentGraph* dg, unsigned person_id, int maternal_allele, int paternal_allele);
    double get_trait_probability(unsigned person_id, enum phased_trait pt);
    
    void preevaluate_init(DescentGraph* dg) {}
    void evaluate_child_peel(unsigned int pmatrix_index, DescentGraph* dg);
    void evaluate_parent_peel(unsigned int pmatrix_index, DescentGraph* dg);
    
 public :
    TraitRfunction(Pedigree* p, GeneticMap* m, unsigned int locus, PeelOperation* po, Rfunction* prev1, Rfunction* prev2) : 
        Rfunction(p, m, locus, po, prev1, prev2),
        half_theta(m->get_theta_halfway(locus)),
        half_inversetheta(m->get_inversetheta_halfway(locus)) {
        
        for(int i = 0; i < 4; ++i) {
            trait_cache[i] = get_trait_probability(peel_id, static_cast<enum phased_trait>(i));
        }
    }
    
    TraitRfunction(const TraitRfunction& rhs) :
        Rfunction(rhs),
        half_theta(rhs.half_theta),
        half_inversetheta(rhs.half_inversetheta) {
    
        for(int i = 0; i < 4; ++i) {
            trait_cache[i] = get_trait_probability(peel_id, static_cast<enum phased_trait>(i));
        }
    }
    
    virtual ~TraitRfunction() {}
    
    TraitRfunction& operator=(const TraitRfunction& rhs) {
        
        if(&rhs != this) {
            Rfunction::operator=(rhs);
            half_theta = rhs.half_theta;
            half_inversetheta = rhs.half_inversetheta;
        }
        
        return *this;
    }
};

#endif

