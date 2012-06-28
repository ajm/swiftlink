#ifndef LKG_SAMPLERRFUNCTION_H_
#define LKG_SAMPLERRFUNCTION_H_

#include <vector>
#include "rfunction.h"

class SamplerRfunction : public Rfunction {
    
    bool ignore_left;
    bool ignore_right;
    vector<double*> transmission;
    vector<unsigned int> children;
    
    
    double get_trait_probability(unsigned person_id, enum phased_trait pt);
    double get_transmission_probability(enum phased_trait parent_trait, enum phased_trait kid_trait, enum parentage parent);
    double get_recombination_probability(DescentGraph* dg, unsigned kid_id, enum phased_trait parent_trait, 
                                         enum phased_trait kid_trait, enum parentage parent);
    void get_recombination_distribution(DescentGraph* dg, unsigned person_id, enum phased_trait parent_trait, 
                                        enum parentage parent, double* dist);
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
    SamplerRfunction(Pedigree* p, GeneticMap* m, unsigned int locus, PeelOperation* po, vector<Rfunction*> previous) : 
        Rfunction(p, m, locus, po, previous), 
        ignore_left(false), 
        ignore_right(false),
        transmission(),
        children() {
        
        set_locus(locus, ignore_left, ignore_right);
        setup_transmission_cache();
        
        for(int i = 0; i < 4; ++i) {
            trait_cache[i] = get_trait_probability(peel_id, static_cast<enum phased_trait>(i));
        }
    }
    
    SamplerRfunction(const SamplerRfunction& rhs) :
        Rfunction(rhs), 
        ignore_left(rhs.ignore_left), 
        ignore_right(rhs.ignore_right),
        transmission(),
        children() {
    
        set_locus(locus, ignore_left, ignore_right);
        setup_transmission_cache();
        
        for(int i = 0; i < 4; ++i) {
            trait_cache[i] = get_trait_probability(peel_id, static_cast<enum phased_trait>(i));
        }
    }
    
    virtual ~SamplerRfunction() {
        teardown_transmission_cache();
    }
    
    SamplerRfunction& operator=(const SamplerRfunction& rhs) {
        
        if(&rhs != this) {
            Rfunction::operator=(rhs);
            ignore_left = rhs.ignore_left;
            ignore_right = rhs.ignore_right;
            
            teardown_transmission_cache();
            setup_transmission_cache();
        }
        
        return *this;
    }
    
    void sample(vector<int>& pmk);
    
    void set_locus(unsigned int l, bool left, bool right) {
        locus = l;
        ignore_left = left;
        ignore_right = right;
        
        theta = theta2 = antitheta = antitheta2 = 1.0;
        
        
        if((locus != 0) and (not ignore_left)) {
            theta2 = map->get_theta(locus-1);
            antitheta2 = map->get_inversetheta(locus-1);
        }
        
        if((locus != (map->num_markers() - 1)) and (not ignore_right)) {
            theta = map->get_theta(locus);
            antitheta = map->get_inversetheta(locus);
        }
        
        for(int i = 0; i < 4; ++i) {
            trait_cache[i] = get_trait_probability(peel_id, static_cast<enum phased_trait>(i));
        }
        
        // XXX this is such a bad idea, it get changed all over the place!
        // but only one thread operates on each locus at a time...
        valid_indices = peel->get_matrix_indices(l);
        
        pmatrix.reset();
        pmatrix_presum.reset();
    }
    
    void set_locus_minimal(unsigned int l) {
        locus = l;
        
        if(locus != 0) {
            theta2     = map->get_theta(locus-1);
            antitheta2 = map->get_inversetheta(locus-1);
        }
        
        if(locus != (map->num_markers() - 1)) {
            theta      = map->get_theta(locus);
            antitheta  = map->get_inversetheta(locus);
        }
        
        for(int i = 0; i < 4; ++i) {
            trait_cache[i] = get_trait_probability(peel_id, static_cast<enum phased_trait>(i));
        }
        
        // XXX this is such a bad idea, it get changed all over the place!
        // but only one thread operates on each locus at a time...
        valid_indices = peel->get_matrix_indices(l);
        
        pmatrix.reset();
        pmatrix_presum.reset();
    }
};

#endif

