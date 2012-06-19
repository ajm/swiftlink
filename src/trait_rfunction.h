#ifndef LKG_TRAITRFUNCTION_H_
#define LKG_TRAITRFUNCTION_H_

#include <vector>
#include "rfunction.h"


class TraitRfunction : public Rfunction {
    
    double get_recombination_probability(DescentGraph* dg, unsigned person_id, int maternal_allele, int paternal_allele);
    double get_trait_probability(unsigned person_id, enum phased_trait pt);
    
    void preevaluate_init(DescentGraph* dg) {}
    void evaluate_child_peel(unsigned int pmatrix_index, DescentGraph* dg);
    void evaluate_parent_peel(unsigned int pmatrix_index, DescentGraph* dg);
    
 public :
    TraitRfunction(Pedigree* p, GeneticMap* m, unsigned int locus, PeelOperation* po, Rfunction* prev1, Rfunction* prev2) : 
        Rfunction(p, m, locus, po, prev1, prev2) {
        
        for(int i = 0; i < 4; ++i) {
            trait_cache[i] = get_trait_probability(peel_id, static_cast<enum phased_trait>(i));
        }
    }
    
    TraitRfunction(const TraitRfunction& rhs) :
        Rfunction(rhs) {
        
        for(int i = 0; i < 4; ++i) {
            trait_cache[i] = get_trait_probability(peel_id, static_cast<enum phased_trait>(i));
        }
    }
    
    virtual ~TraitRfunction() {}
    
    TraitRfunction& operator=(const TraitRfunction& rhs) {
        
        if(&rhs != this) {
            Rfunction::operator=(rhs);
        }
        
        return *this;
    }
    
    void set_thetas(int offset) {
        theta = map->get_theta_partial(locus, offset);
        antitheta = 1.0 - theta;
        
        theta2 = map->get_theta_partial(locus, map->get_lodscore_count() + 1 - offset);
        antitheta2 = 1.0 - theta2;
    }
    
    // for TraitRfunction only trait_cache needs to be set, valid_lod_indices is the
    // same for all loci + pmatrix + pmatrix_presum will not contain different non-zero
    // elements at different loci
    void set_locus_minimal(unsigned int l) {
        locus = l;
        
        for(int i = 0; i < 4; ++i) {
            trait_cache[i] = get_trait_probability(peel_id, static_cast<enum phased_trait>(i));
        }
    }
};

#endif

