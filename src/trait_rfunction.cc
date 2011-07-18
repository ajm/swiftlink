using namespace std;

#include "trait_rfunction.h"
#include "rfunction.h"
#include "descent_graph.h"
#include "genetic_map.h"
#include "pedigree.h"
#include "peeling.h"


TraitRfunction::TraitRfunction(PeelOperation po, Pedigree* p, GeneticMap* m, Rfunction* prev1, Rfunction* prev2) : 
    Rfunction(po, p, m, prev1, prev2) {}

TraitRfunction::TraitRfunction(const TraitRfunction& rhs) :
    Rfunction(rhs) {}
    
TraitRfunction& TraitRfunction::operator=(const TraitRfunction& rhs) {
    
    if(&rhs != this) {
        Rfunction::operator=(rhs);
    }
    
    return *this;
}

// XXX the evaluate method was called with an offset
// and this set as an attribute of the Rfunction super
// class, this needs to be used and passed to the 
// genetic map object
double TraitRfunction::get_recombination_probability(DescentGraph* dg, unsigned locus, unsigned person_id, 
                                                     int maternal_allele, int paternal_allele) {

    double tmp = 1.0;
    double half_recomb_prob;
    
    half_recomb_prob = map->get_theta_halfway(locus);
    
    //printf("HALF-RECOMB: %f\n", half_recomb_prob);
    
    tmp *= ((dg->get(person_id, locus,   MATERNAL) == maternal_allele) ? 1.0 - half_recomb_prob : half_recomb_prob);
    tmp *= ((dg->get(person_id, locus+1, MATERNAL) == maternal_allele) ? 1.0 - half_recomb_prob : half_recomb_prob);
            
    tmp *= ((dg->get(person_id, locus,   PATERNAL) == paternal_allele) ? 1.0 - half_recomb_prob : half_recomb_prob);
    tmp *= ((dg->get(person_id, locus+1, PATERNAL) == paternal_allele) ? 1.0 - half_recomb_prob : half_recomb_prob);
    
    return tmp;
}
    
double TraitRfunction::get_trait_probability(unsigned person_id, enum phased_trait pt, unsigned locus) {
    return (ped->get_by_index(person_id))->get_disease_prob(pt);
}

double TraitRfunction::get_transmission_probability(enum phased_trait parent) {
    return 0.5;
}

