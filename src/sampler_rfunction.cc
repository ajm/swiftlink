using namespace std;

#include <cstdlib>

#include "sampler_rfunction.h"
#include "descent_graph.h"
#include "genetic_map.h"
#include "peeling.h"
#include "peel_matrix.h"


SamplerRfunction::SamplerRfunction(PeelOperation po, Pedigree* p, GeneticMap* m, Rfunction* prev1, Rfunction* prev2) : 
    Rfunction(po, p, m, prev1, prev2) {}

SamplerRfunction::SamplerRfunction(const SamplerRfunction& rhs) :
    Rfunction(rhs) {}

SamplerRfunction& SamplerRfunction::operator=(const SamplerRfunction& rhs) {
    
    if(&rhs != this) {
        Rfunction::operator=(rhs);
    }
    
    return *this;
}

double SamplerRfunction::get_recombination_probability(DescentGraph* dg, unsigned locus, unsigned person_id, 
                                                       int maternal_allele, int paternal_allele) {

    double tmp = 1.0;
    double recomb_prob = 0.0;
    
    if(locus != 0) {
        recomb_prob = exp(map->get_theta(locus-1, temperature));
        tmp *= dg->get(person_id, locus-1, MATERNAL) == maternal_allele ? 1.0 - recomb_prob : recomb_prob;
        tmp *= dg->get(person_id, locus-1, PATERNAL) == paternal_allele ? 1.0 - recomb_prob : recomb_prob;
    }
    
    if(locus != (map->num_markers() - 1)) {
        recomb_prob = exp(map->get_theta(locus, temperature));
        tmp *= dg->get(person_id, locus+1, MATERNAL) == maternal_allele ? 1.0 - recomb_prob : recomb_prob;
        tmp *= dg->get(person_id, locus+1, PATERNAL) == paternal_allele ? 1.0 - recomb_prob : recomb_prob;
    }
    
    return tmp;
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
                //return 0.25;
        }
    //}
    
    return 0.25;
}

double SamplerRfunction::get_transmission_probability(enum phased_trait parent) {
    return ((parent == TRAIT_UU) or (parent == TRAIT_AA)) ? 0.5 : 1.0;
}

double SamplerRfunction::get_transmission_probability2(DescentGraph* dg, 
                                     unsigned locus, 
                                     unsigned person_id, 
                                     enum phased_trait parent_trait, 
                                     enum phased_trait kid_trait, 
                                     enum parentage parent) {
    
    enum trait t = get_trait(kid_trait, parent);
    double tmp = 1.0;
    double recomb_prob;
    
    // deal with homozygotes first
    if(parent_trait == TRAIT_AA) {
        return t == TRAIT_A ? 0.5 : 0.0;
    }
    
    if(parent_trait == TRAIT_UU) {
        return t == TRAIT_U ? 0.5 : 0.0;
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
    
    if(locus != 0) {
        recomb_prob = exp(map->get_theta(locus-1, temperature));
        tmp *= ((dg->get(person_id, locus-1, parent) == p) ? 1.0 - recomb_prob : recomb_prob);
    }
    
    if(locus != (map->num_markers() - 1)) {
        recomb_prob = exp(map->get_theta(locus, temperature));
        tmp *= ((dg->get(person_id, locus+1, parent) == p) ? 1.0 - recomb_prob : recomb_prob);
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
    
    //double recombination_prob;
    double transmission_prob;
    double disease_prob;
    double old_prob1;
    double old_prob2;
    //double tmp = 0.0;
    
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
        
        disease_prob      = get_trait_probability(kid_id, kid_trait, locus);
        //transmission_prob = !dg ?
        //                    0.25 :
        transmission_prob = get_transmission_probability2(dg, locus, kid_id, mat_trait, kid_trait, MATERNAL) *  \
                            get_transmission_probability2(dg, locus, kid_id, pat_trait, kid_trait, PATERNAL);
        
        old_prob1 = previous_rfunction1 != NULL ? previous_rfunction1->get(pmatrix_index) : 1.0;
        old_prob2 = previous_rfunction2 != NULL ? previous_rfunction2->get(pmatrix_index) : 1.0;
        
        pmatrix_presum.set(pmatrix_index, 
                           disease_prob * \
                           transmission_prob * \
                           old_prob1 * \
                           old_prob2);
    }
    
    summation(pmatrix_index, kid_id);
}

