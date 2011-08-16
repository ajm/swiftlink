#ifndef LKG_TRAITRFUNCTION_H_
#define LKG_TRAITRFUNCTION_H_

#include <vector>
#include "rfunction.h"

class TraitRfunction : public Rfunction {
    
    double get_recombination_probability(DescentGraph* dg, unsigned locus, unsigned person_id, 
                                         int maternal_allele, int paternal_allele);
    double get_trait_probability(unsigned person_id, enum phased_trait pt, unsigned locus);
    bool affected_trait(enum phased_trait pt, int allele);
    enum phased_trait get_phased_trait(enum phased_trait m, enum phased_trait p, int maternal_allele, int paternal_allele);
    
    void evaluate_child_peel(PeelMatrixKey& pmatrix_index, DescentGraph* dg, unsigned locus);
    void evaluate_parent_peel(PeelMatrixKey& pmatrix_index, DescentGraph* dg, unsigned locus);
    
 public :
    TraitRfunction(PeelOperation po, Pedigree* p, GeneticMap* m, Rfunction* prev1, Rfunction* prev2);
    virtual ~TraitRfunction() {}
    TraitRfunction(const TraitRfunction& rhs);
    TraitRfunction& operator=(const TraitRfunction& rhs);
    
    double get_result() { 
        return pmatrix.get_result();
    }
};

#endif

