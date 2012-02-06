#ifndef LKG_TRAITRFUNCTION_H_
#define LKG_TRAITRFUNCTION_H_

#include <vector>
#include "rfunction.h"


class TraitRfunction : public Rfunction {
    
    double get_recombination_probability(DescentGraph* dg, unsigned person_id, int maternal_allele, int paternal_allele);
    double get_trait_probability(unsigned person_id, enum phased_trait pt);
    bool affected_trait(enum phased_trait pt, int allele);
    enum phased_trait get_phased_trait(enum phased_trait m, enum phased_trait p, int maternal_allele, int paternal_allele);
    
    void evaluate_child_peel(unsigned int pmatrix_index, DescentGraph* dg);
    void evaluate_parent_peel(unsigned int pmatrix_index, DescentGraph* dg);
    
 public :
    TraitRfunction(Pedigree* p, GeneticMap* m, unsigned int locus, PeelOperation* po, Rfunction* prev1, Rfunction* prev2) : 
        Rfunction(p, m, locus, po, prev1, prev2) {}
    
    TraitRfunction(const TraitRfunction& rhs) :
        Rfunction(rhs) {}
    
    virtual ~TraitRfunction() {}
    
    TraitRfunction& operator=(const TraitRfunction& rhs) {
        
        if(&rhs != this) {
            Rfunction::operator=(rhs);
        }
        
        return *this;
    }
    
    /*
    double get_result() { 
        return pmatrix.get_result();
    }
    */
};

#endif

