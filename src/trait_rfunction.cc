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

double TraitRfunction::get_recombination_probability(DescentGraph* dg, 
                                                     unsigned locus, unsigned person_id, 
                                                     int maternal_allele, int paternal_allele) {
    double tmp = 1.0;    
    double half_theta = map->get_theta_halfway(locus);
    double half_inversetheta = map->get_inversetheta_halfway(locus);
        
    tmp *= ((dg->get(person_id, locus,   MATERNAL) == maternal_allele) ? half_inversetheta : half_theta);
    tmp *= ((dg->get(person_id, locus+1, MATERNAL) == maternal_allele) ? half_inversetheta : half_theta);
            
    tmp *= ((dg->get(person_id, locus,   PATERNAL) == paternal_allele) ? half_inversetheta : half_theta);
    tmp *= ((dg->get(person_id, locus+1, PATERNAL) == paternal_allele) ? half_inversetheta : half_theta);
    
    return tmp;
}
    
double TraitRfunction::get_trait_probability(unsigned person_id, enum phased_trait pt, unsigned locus) {
    return (ped->get_by_index(person_id))->get_disease_prob(pt);
}

bool TraitRfunction::affected_trait(enum phased_trait pt, int allele) {
    
    switch(allele) {
        case 0 :
            return (pt == TRAIT_AU) or (pt == TRAIT_AA);
            
        case 1 :
            return (pt == TRAIT_UA) or (pt == TRAIT_AA);
            
        default :
            abort();
    }
    
    return false;
}

enum phased_trait TraitRfunction::get_phased_trait(enum phased_trait m, enum phased_trait p, 
                                                   int maternal_allele, int paternal_allele) {
    bool m_affected = affected_trait(m, maternal_allele);
    bool p_affected = affected_trait(p, paternal_allele);
    enum phased_trait pt;
    
    if(m_affected) {
        pt = p_affected ? TRAIT_AA : TRAIT_AU;
    }
    else {
        pt = p_affected ? TRAIT_UA : TRAIT_UU;
    }
    
    return pt;
}

/*
double TraitRfunction::get_transmission_probability(enum phased_trait parent) {
    return 0.5;
}
*/

void TraitRfunction::evaluate_child_peel(unsigned int pmatrix_index, 
                                    DescentGraph* dg,
                                    unsigned locus) {
    
    double tmp;
    double total = 0.0;
    
    unsigned int offset = 1 << (2 * peel.get_cutset_size());
    unsigned int presum_index;
    
    enum phased_trait kid_trait;
    enum phased_trait mat_trait;
    enum phased_trait pat_trait;
    
    unsigned kid_id = peel.get_peelnode();
    Person* kid = ped->get_by_index(kid_id);
    
    mat_trait = static_cast<enum phased_trait>((*indices)[pmatrix_index][kid->get_maternalid()]);
    pat_trait = static_cast<enum phased_trait>((*indices)[pmatrix_index][kid->get_paternalid()]);
    
    // iterate over all descent graphs to determine child trait 
    // based on parents' traits
    for(int i = 0; i < 2; ++i) {        // maternal
        for(int j = 0; j < 2; ++j) {    // paternal
            
            kid_trait = get_phased_trait(mat_trait, pat_trait, i, j);
            presum_index = pmatrix_index + (offset * static_cast<int>(kid_trait));
            
            (*indices)[pmatrix_index][kid_id] = static_cast<int>(kid_trait);
            
            tmp = 0.25 * get_trait_probability(kid_id, kid_trait, locus);
            if(tmp == 0.0)
                continue;
            
            tmp *= (!dg ? 1.0 : get_recombination_probability(dg, locus, kid_id, i, j));
            if(tmp == 0.0)
                continue;
            
            tmp *= (previous_rfunction1 != NULL ? previous_rfunction1->get((*indices)[pmatrix_index]) : 1.0) * \
                   (previous_rfunction2 != NULL ? previous_rfunction2->get((*indices)[pmatrix_index]) : 1.0);
            
            pmatrix_presum.add(presum_index, tmp);
            
            total += tmp;
        }
    }
    
    pmatrix.set(pmatrix_index, total);
}

// TODO XXX this is a mess
//
void TraitRfunction::evaluate_parent_peel(
                                     unsigned int pmatrix_index, 
                                     DescentGraph* dg,
                                     unsigned locus) {
    
    unsigned parent_id = peel.get_peelnode();
    unsigned child_node = peel.get_cutnode(0);
    
    unsigned int offset = 1 << (2 * peel.get_cutset_size());
    unsigned int presum_index;
    
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
    enum phased_trait pivot_trait;
    enum phased_trait parent_trait;
    enum phased_trait other_trait;
    
    double tmp;
    double total = 0.0;
    
    other_trait = static_cast<enum phased_trait>((*indices)[pmatrix_index][other_parent_id]);
    
    
    for(unsigned a = 0; a < NUM_ALLELES; ++a) {
        parent_trait = static_cast<enum phased_trait>(a);
        presum_index = pmatrix_index + (offset * a);
        
        (*indices)[pmatrix_index][parent_id] = a;
        
        tmp = get_trait_probability(parent_id, parent_trait, locus);
        if(tmp == 0.0)
            continue;
        
        tmp *= (previous_rfunction1 != NULL ? previous_rfunction1->get((*indices)[pmatrix_index]) : 1.0) * \
               (previous_rfunction2 != NULL ? previous_rfunction2->get((*indices)[pmatrix_index]) : 1.0);
        if(tmp == 0)
            continue;
        
        double child_prob = 1.0;
        
        for(unsigned c = 0; c < peel.get_cutset_size(); ++c) {
            Person* child = ped->get_by_index(peel.get_cutnode(c));
            double child_tmp = 0.0;
            
            if(not child->is_parent(parent_id))
                continue;
            
            for(int i = 0; i < 2; ++i) {        // maternal allele
                for(int j = 0; j < 2; ++j) {    // paternal allele
                    pivot_trait = get_phased_trait(
                                                   ismother ? parent_trait : other_trait, 
                                                   ismother ? other_trait  : parent_trait, 
                                                   i, 
                                                   j);
                    
                    if(pivot_trait != static_cast<enum phased_trait>((*indices)[pmatrix_index][child->get_internalid()]))
                        continue;
                    
                    child_tmp += (!dg ? 0.25 : 0.25 * get_recombination_probability(dg, locus, child->get_internalid(), i, j));
                }
            }
            
            child_prob *= child_tmp;
        }
        
        tmp *= child_prob;
        
        pmatrix_presum.set(presum_index, tmp);
        
        total += tmp;
    }
    
    pmatrix.set(pmatrix_index, total);
}
