#ifndef LKG_SAMPLERRFUNCTION_H_
#define LKG_SAMPLERRFUNCTION_H_

#include <vector>
#include "rfunction.h"

class SamplerRfunction : public Rfunction {
    
    bool ignore_left;
    bool ignore_right;
    
    double get_trait_probability(unsigned person_id, enum phased_trait pt, unsigned locus);
    double get_transmission_probability(enum phased_trait parent_trait, enum phased_trait kid_trait, enum parentage parent);
    double get_recombination_probability(DescentGraph* dg, unsigned locus, unsigned kid_id, 
                                        enum phased_trait parent_trait, enum phased_trait kid_trait, 
                                        enum parentage parent);
                                         
    void evaluate_child_peel(PeelMatrixKey& pmatrix_index, DescentGraph* dg, unsigned locus);
    void evaluate_parent_peel(PeelMatrixKey& pmatrix_index, DescentGraph* dg, unsigned locus);

 public :
    SamplerRfunction(PeelOperation po, Pedigree* p, GeneticMap* m, Rfunction* prev1, Rfunction* prev2);
    virtual ~SamplerRfunction() {}
    SamplerRfunction(const SamplerRfunction& rhs);
    SamplerRfunction& operator=(const SamplerRfunction& rhs);

    void sample(PeelMatrixKey& pmk);
    
    void set_ignore(bool left, bool right) {
        ignore_left = left;
        ignore_right = right;
    }
};

#endif
