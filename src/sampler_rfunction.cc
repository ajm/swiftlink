using namespace std;

#include <cstdlib>

#include "sampler_rfunction.h"
#include "descent_graph.h"
#include "genetic_map.h"
#include "peeling.h"
#include "peel_matrix.h"
#include "random.h"


double SamplerRfunction::get_trait_probability(unsigned person_id, enum phased_trait pt) {
    Person* p = ped->get_by_index(person_id);
    
    if(p->istyped()) {
    
        switch(p->get_marker(locus)) {
            
            case HETERO :
                return ((pt == TRAIT_AU) or (pt == TRAIT_UA)) ? 0.5 : 0.0;
                
            case HOMOZ_A :
                return (pt == TRAIT_UU) ? 1.0 : 0.0;
                
            case HOMOZ_B :
                return (pt == TRAIT_AA) ? 1.0 : 0.0;
                
            default :
                break;
        }
    }
    
    if(not p->isfounder())
              return 0.25;
        
    return map->get_prob(locus, pt);
}

double SamplerRfunction::get_transmission_probability(enum phased_trait parent_trait, enum phased_trait kid_trait, enum parentage parent) {
    enum trait t = get_trait(kid_trait, parent);
    
    if(parent_trait == TRAIT_AA) {
        return (t == TRAIT_A) ? 0.5 : 0.0;
    }
    
    if(parent_trait == TRAIT_UU) {
        return (t == TRAIT_U) ? 0.5 : 0.0;
    }
    
    return 1.0;
}

double SamplerRfunction::get_recombination_probability(DescentGraph* dg, 
                                     unsigned person_id, 
                                     enum phased_trait parent_trait, 
                                     enum phased_trait kid_trait, 
                                     enum parentage parent) {
    
    enum trait t = get_trait(kid_trait, parent);
    
    // deal with homozygotes first
    if(parent_trait == TRAIT_AA) {
        //return (t == TRAIT_A) ? 0.5 : 0.0;
        return (t == TRAIT_A) ? 0.25 : 0.0; // XXX <--- transmission prob as well!!!
    }
    else if(parent_trait == TRAIT_UU) {
        //return (t == TRAIT_U) ? 0.5 : 0.0;
        return (t == TRAIT_U) ? 0.25 : 0.0; // XXX <--- transmission prob as well!!!
    }
    
    // heterozygotes are informative, so i can look up
    // the recombination fractions
    int p = 0;
    if(parent_trait == TRAIT_UA) {
        p = (t == TRAIT_U) ? 0 : 1;
    }
    else if(parent_trait == TRAIT_AU) {
        p = (t == TRAIT_A) ? 0 : 1;
    }
    
    double tmp = 0.5; // XXX <--- transmission prob
    if((locus != 0) and (not ignore_left)) {
        tmp *= ((dg->get(person_id, locus-1, parent) == p) ? antitheta2 : theta2);
    }
    
    if((locus != (map->num_markers() - 1)) and (not ignore_right)) {
        tmp *= ((dg->get(person_id, locus+1, parent) == p) ? antitheta : theta);
    }
    
    return tmp;
}

void SamplerRfunction::sample(vector<int>& pmk) {
    double prob_dist[NUM_ALLELES];
    double total = 0.0;
    unsigned node = peel->get_peelnode();
    enum phased_trait trait;
    unsigned int last = 0;
    
    // extract probabilities
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        trait = static_cast<enum phased_trait>(i);
        pmk[node] = i;
        
        prob_dist[i] = pmatrix_presum.get(pmk);
        total += prob_dist[i];        
    }
    
    if(total == 0.0) {
        fprintf(stderr, "error: probabilities sum to zero\n");
        abort();
    }
    
    // normalise
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        prob_dist[i] /= total;
    }
    
    // sample
    double r = get_random();
    total = 0.0;
    
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        total += prob_dist[i];
        if(r <= total) {
            pmk[node] = i;
            return;
        }
        
        if(prob_dist[i] != 0.0) {
            last = i;
        }
    }
    
    // XXX still abort?
    abort();
}

void SamplerRfunction::evaluate_child_peel(unsigned int pmatrix_index, DescentGraph* dg) {
        
    unsigned int presum_index;
    Person* kid = ped->get_by_index(peel_id);    
    
    enum phased_trait mat_trait;
    enum phased_trait pat_trait;
    enum phased_trait kid_trait;
    
    double tmp;
    double total = 0.0;
    
    
    mat_trait = static_cast<enum phased_trait>(indices[pmatrix_index][kid->get_maternalid()]);
    pat_trait = static_cast<enum phased_trait>(indices[pmatrix_index][kid->get_paternalid()]);
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        kid_trait = static_cast<enum phased_trait>(i);
        presum_index = pmatrix_index + (index_offset * i);
        
        indices[pmatrix_index][peel_id] = i;
        
        tmp = get_trait_probability(peel_id, kid_trait);
        if(tmp == 0.0)
            continue;
        
        tmp *= (get_recombination_probability(dg, peel_id, mat_trait, kid_trait, MATERNAL) * \
                get_recombination_probability(dg, peel_id, pat_trait, kid_trait, PATERNAL));
        if(tmp == 0.0)
            continue;
        
        tmp *= ((previous_rfunction1 != NULL) ? previous_rfunction1->get(indices[pmatrix_index]) : 1.0) * \
               ((previous_rfunction2 != NULL) ? previous_rfunction2->get(indices[pmatrix_index]) : 1.0);
        
        pmatrix_presum.set(presum_index, tmp);
        
        total += tmp;
    }
    
    pmatrix.set(pmatrix_index, total);
}

void SamplerRfunction::evaluate_parent_peel(unsigned int pmatrix_index, DescentGraph* dg) {
    
    unsigned int presum_index;
    
    enum phased_trait child_trait;
    enum phased_trait parent_trait;
    enum phased_trait other_trait;
    
    double tmp;
    double total = 0.0;
    
    
    for(unsigned int i = 0; i < NUM_ALLELES; ++i) {
        parent_trait = static_cast<enum phased_trait>(i);
        presum_index = pmatrix_index + (index_offset * i);
        
        indices[pmatrix_index][peel_id] = i;
        
        tmp = get_trait_probability(peel_id, parent_trait);
        if(tmp == 0.0)
            continue;
        
        tmp *= ((previous_rfunction1 != NULL ? previous_rfunction1->get(indices[pmatrix_index]) : 1.0) * \
                (previous_rfunction2 != NULL ? previous_rfunction2->get(indices[pmatrix_index]) : 1.0));
        if(tmp == 0.0)
            continue;
        
        double child_prob = 1.0;
        
        for(unsigned int c = 0; c < peel->get_cutset_size(); ++c) {
            unsigned int child_id = peel->get_cutnode(c);
            Person* child = ped->get_by_index(child_id);
            
            if(not child->is_parent(peel_id))
                continue;
            
            child_trait = static_cast<enum phased_trait>(indices[pmatrix_index][child_id]);
            other_trait = static_cast<enum phased_trait>(indices[pmatrix_index][child->get_maternalid() == peel_id ? child->get_paternalid() : child->get_maternalid()]);
            
            if(child->get_maternalid() == peel_id) {
                child_prob *= (get_recombination_probability(dg, child_id, parent_trait, child_trait, MATERNAL) *  \
                               get_recombination_probability(dg, child_id, other_trait, child_trait, PATERNAL));
            }
            else {
                child_prob *= (get_recombination_probability(dg, child_id, parent_trait, child_trait, PATERNAL) *  \
                               get_recombination_probability(dg, child_id, other_trait,  child_trait, MATERNAL));
            }
            
            if(child_prob == 0.0)
                break;
        }
        
        tmp *= child_prob;
        
        pmatrix_presum.set(presum_index, tmp);
        
        total += tmp;
    }
    
    pmatrix.set(pmatrix_index, total);
}
