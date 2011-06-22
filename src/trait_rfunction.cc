using namespace std;

#include <vector>

#include "trait_rfunction.h"
#include "rfunction.h"
#include "descent_graph.h"
#include "genetic_map.h"
#include "pedigree.h"


TraitRfunction::TraitRfunction(PeelOperation po, Pedigree* p, GeneticMap* m, 
                               vector<Rfunction*>& previous_functions, unsigned index) : 
    Rfunction(po, p, m, previous_functions, index) {}

TraitRfunction::TraitRfunction(const Rfunction& rhs) :
    Rfunction(rhs) {}
    
TraitRfunction& TraitRfunction::operator=(const Rfunction& rhs) {
    
    if(&rhs != this) {
        Rfunction::operator=(rhs);
    }
    
    return *this;
}

// XXX this needs to happen at arbitrary, user-defined (?) increments
// along the space between two markers
// i need to keep this interface the same to benefit from inheritance
// so the offset will need to be set elsewhere
// --> evaluate method has it as a parameter, so it should be easy enough... :-P
double TraitRfunction::get_recombination_probability(DescentGraph* dg, unsigned locus, unsigned person_id, 
                                                     int maternal_allele, int paternal_allele) {

    double tmp = 1.0;
    double half_recomb_prob;
    
    half_recomb_prob = map->get_theta_halfway(locus);
    
    tmp *= dg->get(person_id, locus,   MATERNAL) == maternal_allele ? 1.0 - half_recomb_prob : half_recomb_prob;
    tmp *= dg->get(person_id, locus+1, MATERNAL) == maternal_allele ? 1.0 - half_recomb_prob : half_recomb_prob;
            
    tmp *= dg->get(person_id, locus,   PATERNAL) == paternal_allele ? 1.0 - half_recomb_prob : half_recomb_prob;
    tmp *= dg->get(person_id, locus+1, PATERNAL) == paternal_allele ? 1.0 - half_recomb_prob : half_recomb_prob;
    
    return tmp;
}
    
double TraitRfunction::get_trait_probability(unsigned person_id, enum phased_trait pt, unsigned locus) {
    return (ped->get_by_index(person_id))->get_disease_prob(pt);
}

