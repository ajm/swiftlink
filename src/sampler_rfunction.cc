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
    
    // penetrace prob
    if(not p->isfounder()) {
        if(p->istyped()) {
            switch(p->get_marker(locus)) {
                case HETERO :
                    return ((pt == TRAIT_AU) or (pt == TRAIT_UA)) ? 1.0 : 0.0;
                case HOMOZ_A :
                    return (pt == TRAIT_UU) ? 1.0 : 0.0;
                case HOMOZ_B :
                    return (pt == TRAIT_AA) ? 1.0 : 0.0;
                default :
                    return 1.0;
            }
        }
        else {
            return 1.0;
        }
    }
    
    // penetrance + founder prob
    if(p->istyped()) {
        switch(p->get_marker(locus)) {
            case HETERO :
                return ((pt == TRAIT_AU) or (pt == TRAIT_UA)) ? map->get_prob(locus, pt) : 0.0;
            case HOMOZ_A :
                return (pt == TRAIT_UU) ? map->get_prob(locus, pt) : 0.0;
            case HOMOZ_B :
                return (pt == TRAIT_AA) ? map->get_prob(locus, pt) : 0.0;
            default :
                return map->get_prob(locus, pt);
        }
    }
    else {
        return map->get_prob(locus, pt);
    }
    
    fprintf(stderr, "error: %s:%d\n", __FILE__, __LINE__);
    abort();
}

double SamplerRfunction::get_recombination_probability(DescentGraph* dg, 
                                     unsigned person_id, 
                                     enum phased_trait parent_trait, 
                                     enum phased_trait kid_trait, 
                                     enum parentage parent) {
    
    enum trait t = get_trait(kid_trait, parent);
    
    // recombination + transmission prob
    double tmp = 1.0;
    
    // deal with homozygotes first
    if(parent_trait == TRAIT_AA) {
        tmp = (t == TRAIT_A) ? 1.0 : 0.0;
        if((locus != 0) and (not ignore_left))
            tmp *= 0.5;
        if((locus != (map->num_markers() - 1)) and (not ignore_right))
            tmp *= 0.5;
        return tmp;
    }
    else if(parent_trait == TRAIT_UU) {
        tmp = (t == TRAIT_U) ? 1.0 : 0.0;
        if((locus != 0) and (not ignore_left))
            tmp *= 0.5;
        if((locus != (map->num_markers() - 1)) and (not ignore_right))
            tmp *= 0.5;
        return tmp;
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
    
    tmp = 0.5;
    
    if((locus != 0) and (not ignore_left)) {
        tmp *= ((dg->get(person_id, locus-1, parent) == p) ? antitheta2 : theta2);
    }
    
    if((locus != (map->num_markers() - 1)) and (not ignore_right)) {
        tmp *= ((dg->get(person_id, locus+1, parent) == p) ? antitheta : theta);
    }
    
    return tmp;
}

void SamplerRfunction::sample(DescentGraph* dg, vector<int>& pmk) {
    double prob_dist[NUM_ALLELES];
    
    /*
    // get probabilities
    switch(peel->get_type()) {
        case CHILD_PEEL:
            sample_child(dg, pmk, prob_dist);
            break;
        case LAST_PEEL:
        case PARTNER_PEEL:
            sample_partner(pmk, prob_dist);
            break;
        default:
            fprintf(stderr, "i hav' not written the code for this!\n");
            abort();
    }
    */
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        pmk[peel_id] = i;
        prob_dist[i] = pmatrix_presum.get(pmk);        
    }
    
    normalise(prob_dist, 4);
    
        
    // sample
    unsigned int last = 0;
    double r = get_random();
    double total = 0.0;
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        total += prob_dist[i];
        if(r < total) {
            pmk[peel_id] = i;
            return;
        }
        
        if(prob_dist[i] != 0.0) {
            last = i;
        }
    }
    
    pmk[peel_id] = last;
}

void SamplerRfunction::sample_child(DescentGraph* dg, vector<int>& pmk, double* prob) {
    Person* kid = ped->get_by_index(peel_id);    
    
    enum phased_trait mat_trait;
    enum phased_trait pat_trait;
    enum phased_trait kid_trait;
    
    double tmp;
    double trans[4];
    double prev1[4];
    double prev2[4];
    double total;
    
    mat_trait = static_cast<enum phased_trait>(pmk[kid->get_maternalid()]);
    pat_trait = static_cast<enum phased_trait>(pmk[kid->get_paternalid()]);
    
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        kid_trait = static_cast<enum phased_trait>(i);
        
        trans[i] = (get_recombination_probability(dg, peel_id, mat_trait, kid_trait, MATERNAL) * \
                    get_recombination_probability(dg, peel_id, pat_trait, kid_trait, PATERNAL));
    }
    
    normalise(trans, 4);
    
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        pmk[peel_id] = i;
        
        prev1[i] = ((previous_rfunction1 != NULL) ? previous_rfunction1->get(pmk) : 1.0);
        prev2[i] = ((previous_rfunction2 != NULL) ? previous_rfunction2->get(pmk) : 1.0);
    }
    
    //normalise(prev1, 4);
    //normalise(prev2, 4);

    total = 0.0;
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        kid_trait = static_cast<enum phased_trait>(i);

        tmp = get_trait_probability(peel_id, kid_trait) * trans[i] * prev1[i] * prev2[i];

        prob[i] = tmp;
        total += tmp;
    }
    
    if(total == 0.0) {
        fprintf(stderr, "error: probabilities sum to zero (%s:%d)\n", __FILE__, __LINE__);
        abort();
    }
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        prob[i] /= total;
    }
}

void SamplerRfunction::sample_partner(vector<int>& pmk, double* prob) {
    enum phased_trait partner_trait;    
    double tmp;
    double total;    
    double prev1[4];
    double prev2[4];
    
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        pmk[peel_id] = i;
        prev1[i] = ((previous_rfunction1 != NULL) ? previous_rfunction1->get(pmk) : 1.0);
        prev2[i] = ((previous_rfunction2 != NULL) ? previous_rfunction2->get(pmk) : 1.0);
    }
    
    //normalise(prev1, 4);
    //normalise(prev2, 4);
    
    total = 0.0;
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        partner_trait = static_cast<enum phased_trait>(i);
        
        tmp = get_trait_probability(peel_id, partner_trait) * prev1[i] * prev2[i];
        
        prob[i] = tmp;
        total += tmp;
    }
    
    if(total == 0.0) {
        fprintf(stderr, "error: probabilities sum to zero (%d)\n", peel_id);
        abort();
    }
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        prob[i] /= total;
    }
}

void SamplerRfunction::evaluate_child_peel(unsigned int pmatrix_index, DescentGraph* dg) {
        
    unsigned int presum_index;
    Person* kid = ped->get_by_index(peel_id);    
    
    enum phased_trait mat_trait;
    enum phased_trait pat_trait;
    enum phased_trait kid_trait;
    double tmp;
    double total = 0.0;
    
    double trans[4];
    double pen[4];
    
    mat_trait = static_cast<enum phased_trait>(indices[pmatrix_index][kid->get_maternalid()]);
    pat_trait = static_cast<enum phased_trait>(indices[pmatrix_index][kid->get_paternalid()]);
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        kid_trait = static_cast<enum phased_trait>(i);
        trans[i] = get_recombination_probability(dg, peel_id, mat_trait, kid_trait, MATERNAL) * \
                   get_recombination_probability(dg, peel_id, pat_trait, kid_trait, PATERNAL);
        pen[i] = get_trait_probability(peel_id, kid_trait);
    }
    
    normalise(trans, 4);
    //normalise(pen, 4);
    
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        kid_trait = static_cast<enum phased_trait>(i);
        presum_index = pmatrix_index + (index_offset * i);
        
        indices[pmatrix_index][peel_id] = i;
        
        //tmp = get_trait_probability(peel_id, kid_trait);
        tmp = pen[i];
        if(tmp == 0.0)
            continue;
        
        //tmp *= (get_recombination_probability(dg, peel_id, mat_trait, kid_trait, MATERNAL) *
        //        get_recombination_probability(dg, peel_id, pat_trait, kid_trait, PATERNAL));
        tmp *= trans[i];
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
    
    enum phased_trait pivot_trait;
    enum phased_trait mat_trait;
    enum phased_trait pat_trait;
    
    double tmp;
    double total = 0.0;
    
    
    for(unsigned int i = 0; i < NUM_ALLELES; ++i) {
        mat_trait = pat_trait = static_cast<enum phased_trait>(i);
        presum_index = pmatrix_index + (index_offset * i);
        
        indices[pmatrix_index][peel_id] = i;
        
        tmp = get_trait_probability(peel_id, mat_trait);
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
            
            
            pivot_trait = static_cast<enum phased_trait>(indices[pmatrix_index][child_id]);
            
            if(child->get_maternalid() == peel_id) {
                pat_trait = static_cast<enum phased_trait>(indices[pmatrix_index][child->get_paternalid()]);
            }
            else {
                mat_trait = static_cast<enum phased_trait>(indices[pmatrix_index][child->get_maternalid()]);
            }
            
            
            child_prob *= (get_recombination_probability(dg, child_id, mat_trait, pivot_trait, MATERNAL) *  \
                           get_recombination_probability(dg, child_id, pat_trait, pivot_trait, PATERNAL));
            
            if(child_prob == 0.0)
                break;
        }
        
        tmp *= child_prob;
        
        pmatrix_presum.set(presum_index, tmp);
        
        total += tmp;
    }
    
    pmatrix.set(pmatrix_index, total);
}

