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
    double homoz_cache;
    vector<double*> transmission;
    vector<double> thetas;
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
    SamplerRfunction(Pedigree* p, GeneticMap* m, unsigned int locus, PeelOperation* po, Rfunction* prev1, Rfunction* prev2) : 
        Rfunction(p, m, locus, po, prev1, prev2), 
        ignore_left(false), 
        ignore_right(false),
        theta(0.0),
        antitheta(1.0),
        theta2(0.0),
        antitheta2(1.0),
        homoz_cache(1.0),
        transmission(),
        thetas(4, 0.0),
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
        theta(rhs.theta),
        antitheta(rhs.antitheta),
        theta2(rhs.theta2),
        antitheta2(rhs.antitheta2),
        homoz_cache(rhs.homoz_cache),
        transmission(),
        thetas(4, 0.0),
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
            theta = rhs.theta;
            antitheta = rhs.antitheta;
            theta2 = rhs.theta2;
            antitheta2 = rhs.antitheta2;
            homoz_cache = rhs.homoz_cache;
            
            thetas = rhs.thetas;
            
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
        
        homoz_cache = 1.0;
        
        
        if((locus != 0) and (not ignore_left)) {
            theta2 = map->get_theta(locus-1);
            antitheta2 = map->get_inversetheta(locus-1);
            homoz_cache *= 0.5;
        }
        
        if((locus != (map->num_markers() - 1)) and (not ignore_right)) {
            theta = map->get_theta(locus);
            antitheta = map->get_inversetheta(locus);
            homoz_cache *= 0.5;
        }
        
        for(int i = 0; i < 4; ++i) {
            trait_cache[i] = get_trait_probability(peel_id, static_cast<enum phased_trait>(i));
        }
        
        // XXX this is such a bad idea, it get changed all over the place!
        valid_indices = peel->get_valid_indices(l);
        
        pmatrix.reset();
        pmatrix_presum.reset();
    }
};

#endif

