#ifndef LKG_SAMPLERRFUNCTION_H_
#define LKG_SAMPLERRFUNCTION_H_

#include <vector>
#include "rfunction.h"

class SamplerRfunction : public Rfunction {
    
    bool ignore_left;
    bool ignore_right;
    double theta;
    double antitheta;
    double theta2;
    double antitheta2;
    vector<double*> transmission;
    vector<unsigned int> children;
    
    
    double get_trait_probability(unsigned person_id, enum phased_trait pt);
    double get_transmission_probability(enum phased_trait parent_trait, enum phased_trait kid_trait, enum parentage parent);
    double get_recombination_probability(DescentGraph* dg, unsigned kid_id, enum phased_trait parent_trait, 
                                         enum phased_trait kid_trait, enum parentage parent);
    void transmission_matrix(DescentGraph* dg, int kid_id, double* tmatrix);
    void populate_transmission_cache(DescentGraph* dg);
    void setup_transmission_cache();
    void teardown_transmission_cache();
    unsigned int transmission_index(enum phased_trait mat_trait, 
                                    enum phased_trait pat_trait, 
                                    enum phased_trait kid_trait);
    void evaluate_child_peel(unsigned int pmatrix_index, DescentGraph* dg);
    void evaluate_parent_peel(unsigned int pmatrix_index, DescentGraph* dg);
    
    void preevaluate_init(DescentGraph* dg);

 public :    
    SamplerRfunction(Pedigree* p, GeneticMap* m, unsigned int locus, PeelOperation* po, Rfunction* prev1, Rfunction* prev2) : 
        Rfunction(p, m, locus, po, prev1, prev2), 
        ignore_left(false), 
        ignore_right(false),
        theta(0.0),
        antitheta(0.0),
        theta2(0.0),
        antitheta2(0.0),
        transmission(),
        children() {
        
        set_locus(locus);
        setup_transmission_cache();
    }
    
    SamplerRfunction(const SamplerRfunction& rhs) :
        Rfunction(rhs), 
        ignore_left(rhs.ignore_left), 
        ignore_right(rhs.ignore_right),
        theta(rhs.theta),
        antitheta(rhs.antitheta),
        theta2(rhs.theta2),
        antitheta2(rhs.antitheta2),
        transmission(),
        children() {
    
        set_locus(locus);
        setup_transmission_cache();    
    }
    
    virtual ~SamplerRfunction() {
        teardown_transmission_cache();
    }
    
    SamplerRfunction& operator=(const SamplerRfunction& rhs) {
        
        if(&rhs != this) {
            Rfunction::operator=(rhs);
            ignore_left = rhs.ignore_left;
            ignore_right = rhs.ignore_right;
            theta = rhs.theta;
            antitheta = rhs.antitheta;
            theta2 = rhs.theta2;
            antitheta2 = rhs.antitheta2;
            
            teardown_transmission_cache();
            setup_transmission_cache();
        }
        
        return *this;
    }
    
    void sample(vector<int>& pmk);
    
    void set_ignores(bool left, bool right) {
        ignore_left = left;
        ignore_right = right;
    }
    
    void set_locus(unsigned int l) {
        locus = l;
        
        if(locus != (map->num_markers() - 1)) {
            theta = map->get_theta(locus);
            antitheta = map->get_inversetheta(locus);
        }
        
        if(locus != 0) {
            theta2 = map->get_theta(locus-1);
            antitheta2 = map->get_inversetheta(locus-1);
        }
    }
};

#endif

