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
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        total += prob_dist[i];
        if(r < total) {
            pmk.add(node, static_cast<enum phased_trait>(i));
            return;
        }
    }
    
    abort();
}

