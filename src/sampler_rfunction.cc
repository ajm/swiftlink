using namespace std;

#include <cstdlib>

#include "sampler_rfunction.h"
#include "descent_graph.h"
#include "genetic_map.h"
#include "peeling.h"
#include "peel_matrix.h"


SamplerRfunction::SamplerRfunction(PeelOperation po, Pedigree* p, GeneticMap* m, Rfunction* prev1, Rfunction* prev2) : 
    Rfunction(po, p, m, prev1, prev2), 
    ignore_left(false), 
    ignore_right(false) {}

SamplerRfunction::SamplerRfunction(const SamplerRfunction& rhs) :
    Rfunction(rhs), 
    ignore_left(rhs.ignore_left), 
    ignore_right(rhs.ignore_right) {}

SamplerRfunction& SamplerRfunction::operator=(const SamplerRfunction& rhs) {
    
    if(&rhs != this) {
        Rfunction::operator=(rhs);
        ignore_left = rhs.ignore_left;
        ignore_right = rhs.ignore_right;
    }
    
    return *this;
}
    
double SamplerRfunction::get_trait_probability(unsigned person_id, enum phased_trait pt, unsigned locus) {
    Person* p = ped->get_by_index(person_id);
    
    //if(p->istyped()) {
    
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
    //}
    
    // these need to be properly normalised for each marker in the map, right?
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
                                     unsigned locus, 
                                     unsigned person_id, 
                                     enum phased_trait parent_trait, 
                                     enum phased_trait kid_trait, 
                                     enum parentage parent) {
    
    enum trait t = get_trait(kid_trait, parent);
    double tmp = 1.0;
    
    // deal with homozygotes first
    if(parent_trait == TRAIT_AA) {
        return (t == TRAIT_A) ? 0.5 : 0.0;
    }
    
    if(parent_trait == TRAIT_UU) {
        return (t == TRAIT_U) ? 0.5 : 0.0;
    }
    
    // heterozygotes are informative, so i can look up
    // the recombination fractions
    char p = 0;
    if(parent_trait == TRAIT_UA) {
        p = (t == TRAIT_U) ? 0 : 1;
    }
    else if(parent_trait == TRAIT_AU) {
        p = (t == TRAIT_A) ? 0 : 1;
    }
    
    if((locus != 0) and (not ignore_left)) {
        tmp *= ((dg->get(person_id, locus-1, parent) == p) ? map->get_inversetheta(locus-1) : map->get_theta(locus-1));
    }
    
    if((locus != (map->num_markers() - 1)) and (not ignore_right)) {
        tmp *= ((dg->get(person_id, locus+1, parent) == p) ? map->get_inversetheta(locus) : map->get_theta(locus));
    }
    
    return tmp;
}

void SamplerRfunction::sample(PeelMatrixKey& pmk) {
    double prob_dist[NUM_ALLELES];
    double total = 0.0;
    unsigned node = peel.get_peelnode();
    enum phased_trait trait;
    
    // extract probabilities
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        trait = static_cast<enum phased_trait>(i);
        pmk.add(node, trait);
        
        prob_dist[i] = pmatrix_presum.get(pmk);
        total += prob_dist[i];        
    }
    
    // normalise
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        prob_dist[i] /= total;
    }
    
    // sample
    double r = random() / static_cast<double>(RAND_MAX);
    total = 0.0;
    
    /*
    printf("random = %f\n", r);
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        printf("  sample[%d] = %f\n", i, prob_dist[i]);
    }
    */
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        total += prob_dist[i];
        if(r < total) {
            pmk.add(node, static_cast<enum phased_trait>(i));
            return;
        }
    }
    
    abort();
}

void SamplerRfunction::evaluate_child_peel(
                    PeelMatrixKey& pmatrix_index, 
                    DescentGraph* dg,
                    unsigned locus) {
    
    double transmission_prob;
    double recombination_prob;
    double disease_prob;
    double old_prob1;
    double old_prob2;
    
    enum phased_trait kid_trait;
    enum phased_trait mat_trait;
    enum phased_trait pat_trait;
    
    unsigned kid_id = peel.get_peelnode();
    Person* kid = ped->get_by_index(kid_id);
    
    mat_trait = pmatrix_index.get(kid->get_maternalid());
    pat_trait = pmatrix_index.get(kid->get_paternalid());
    
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        kid_trait = static_cast<enum phased_trait>(i);
        pmatrix_index.add(kid_id, kid_trait);
        
        disease_prob       = get_trait_probability(kid_id, kid_trait, locus);
        transmission_prob  = get_transmission_probability(mat_trait, kid_trait, MATERNAL) * \
                             get_transmission_probability(pat_trait, kid_trait, PATERNAL);
        recombination_prob = get_recombination_probability(dg, locus, kid_id, mat_trait, kid_trait, MATERNAL) *  \
                             get_recombination_probability(dg, locus, kid_id, pat_trait, kid_trait, PATERNAL);
        
        old_prob1 = previous_rfunction1 != NULL ? previous_rfunction1->get(pmatrix_index) : 1.0;
        old_prob2 = previous_rfunction2 != NULL ? previous_rfunction2->get(pmatrix_index) : 1.0;
        
        pmatrix_presum.set(pmatrix_index, \
                           disease_prob * \
                           transmission_prob * \
                           recombination_prob * \
                           old_prob1 * \
                           old_prob2);
    }
    
    summation(pmatrix_index, kid_id);
}

// TODO XXX this is a mess
//
void SamplerRfunction::evaluate_parent_peel(
                                          PeelMatrixKey& pmatrix_index, 
                                          DescentGraph* dg,
                                          unsigned locus) {
    double disease_prob;
    double old_prob1;
    double old_prob2;    
    
    unsigned parent_id = peel.get_peelnode();
    unsigned child_node = peel.get_cutnode(0);
    
    // find a child of parent_id
    for(unsigned i = 0; i < peel.get_cutset_size(); ++i) {
        Person* ptmp = ped->get_by_index(peel.get_cutnode(i));
        if(ptmp->is_parent(parent_id)) {
            child_node = peel.get_cutnode(i);
            break;
        }
    }
    
    Person* p = ped->get_by_index(child_node);
    bool ismother = parent_id == p->get_maternalid();
    unsigned other_parent_id = ismother ? p->get_paternalid() : p->get_maternalid();
    enum phased_trait child_trait;
    enum phased_trait parent_trait;
    enum phased_trait other_trait;
    double tmp = 0.0;
    
    other_trait = pmatrix_index.get(other_parent_id);
    
    
    for(unsigned a = 0; a < NUM_ALLELES; ++a) {
        parent_trait = static_cast<enum phased_trait>(a);
        pmatrix_index.add(parent_id, parent_trait);
        
        disease_prob = get_trait_probability(parent_id, parent_trait, locus);
        
        old_prob1 = previous_rfunction1 != NULL ? previous_rfunction1->get(pmatrix_index) : 1.0;
        old_prob2 = previous_rfunction2 != NULL ? previous_rfunction2->get(pmatrix_index) : 1.0;
        
        double child_prob = 1.0;
        
        for(unsigned c = 0; c < peel.get_cutset_size(); ++c) {
            Person* child = ped->get_by_index(peel.get_cutnode(c));
            
            if(not child->is_parent(parent_id))
                continue;
            
            child_trait = pmatrix_index.get(child->get_internalid());
            
            child_prob *=  (get_transmission_probability(parent_trait, child_trait, ismother ? MATERNAL : PATERNAL) * \
                            get_transmission_probability(other_trait,  child_trait, ismother ? PATERNAL : MATERNAL) * \
                            get_recombination_probability(dg, locus, child->get_internalid(), parent_trait, child_trait, ismother ? MATERNAL : PATERNAL) *  \
                            get_recombination_probability(dg, locus, child->get_internalid(), other_trait,  child_trait, ismother ? PATERNAL : MATERNAL));            
        }
        
        tmp = (child_prob * disease_prob * old_prob1 * old_prob2);
        
        pmatrix_presum.set(pmatrix_index, tmp);
    }
    
    summation(pmatrix_index, parent_id);
}
