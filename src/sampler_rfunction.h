#ifndef LKG_SAMPLERRFUNCTION_H_
#define LKG_SAMPLERRFUNCTION_H_

#include <vector>
#include "rfunction.h"

class SamplerRfunction : public Rfunction {
    
    bool ignore_left;
    bool ignore_right;
    
    double get_trait_probability(unsigned person_id, enum phased_trait pt);
    double get_transmission_probability(enum phased_trait parent_trait, enum phased_trait kid_trait, enum parentage parent);
    double get_recombination_probability(DescentGraph* dg, unsigned kid_id, enum phased_trait parent_trait, 
                                         enum phased_trait kid_trait, enum parentage parent);
    
    void evaluate_child_peel(unsigned int pmatrix_index, DescentGraph* dg);
    void evaluate_parent_peel(unsigned int pmatrix_index, DescentGraph* dg);
    //void evaluate_partner_peel(unsigned int pmatrix_index);
    
    void sample_child(DescentGraph* dg, vector<int>& pmk, double* prob);
    void sample_partner(vector<int>& pmk, double* prob);
    

 public :    
    SamplerRfunction(Pedigree* p, GeneticMap* m, unsigned int locus, PeelOperation* po, Rfunction* prev1, Rfunction* prev2) : 
        Rfunction(p, m, locus, po, prev1, prev2), 
        ignore_left(false), 
        ignore_right(false) {}
    
    SamplerRfunction(const SamplerRfunction& rhs) :
        Rfunction(rhs), 
        ignore_left(rhs.ignore_left), 
        ignore_right(rhs.ignore_right) {}
    
    virtual ~SamplerRfunction() {}
    
    SamplerRfunction& operator=(const SamplerRfunction& rhs) {
        
        if(&rhs != this) {
            Rfunction::operator=(rhs);
            ignore_left = rhs.ignore_left;
            ignore_right = rhs.ignore_right;
        }
        
        return *this;
    }
    
    void sample(DescentGraph* dg, vector<int>& pmk);
    
    void set_ignores(bool left, bool right) {
        ignore_left = left;
        ignore_right = right;
    }
};

#endif

