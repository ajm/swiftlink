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

/*
double SamplerRfunction::get_recombination_probability(DescentGraph* dg, 
                                     unsigned person_id, 
                                     enum phased_trait parent_trait, 
                                     enum phased_trait kid_trait, 
                                     enum parentage parent) {
    
    enum trait t = get_trait(kid_trait, parent);
    int p = 0;
    double tmp = 1.0;
    
    // recombination + transmission prob
    
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
    if(parent_trait == TRAIT_UA) {
        p = (t == TRAIT_U) ? 0 : 1;
    }
    else if(parent_trait == TRAIT_AU) {
        p = (t == TRAIT_A) ? 0 : 1;
    }
    
    tmp = 0.5;
    
    
    //bool crossover1 = false;
    //bool crossover2 = false;
        
    if((locus != 0) and (not ignore_left)) {
        tmp *= ((dg->get(person_id, locus-1, parent) == p) ? antitheta2 : theta2);
        //crossover1 = dg->get(person_id, locus-1, parent) != p;
    }
    
    if((locus != (map->num_markers() - 1)) and (not ignore_right)) {
        tmp *= ((dg->get(person_id, locus+1, parent) == p) ? antitheta : theta);
        //crossover2 = dg->get(person_id, locus+1, parent) != p;
    }
    
    
    //if(crossover1) {
    //    tmp *= (crossover2 ? thetas[0] : thetas[1]);
    //}
    //else {
    //    tmp *= (crossover2 ? thetas[2] : thetas[3]);
    //}
    
    
    return tmp;
}
*/

double SamplerRfunction::get_recombination_probability(DescentGraph* dg, 
                                     unsigned person_id, 
                                     enum phased_trait parent_trait, 
                                     enum phased_trait kid_trait, 
                                     enum parentage parent) {
    
    enum trait t = get_trait(kid_trait, parent);
    int p = 0;
    
    switch(parent_trait) {
        case TRAIT_UU :
            return (t == TRAIT_U) ? 1.0 /*homoz_cache*/ : 0.0;
    
        case TRAIT_AU :
            p = (t == TRAIT_A) ? 0 : 1;
            break;
        
        case TRAIT_UA :
            p = (t == TRAIT_U) ? 0 : 1;
            break;
        
        case TRAIT_AA :
            return (t == TRAIT_A) ? 1.0 /*homoz_cache*/ : 0.0;
    }
    
    double tmp = 0.5;
    
    if(locus != 0) {
        tmp *= ((dg->get(person_id, locus-1, parent) == p) ? antitheta2 : theta2);
    }
    
    if(locus != (map->num_markers() - 1)) {
        tmp *= ((dg->get(person_id, locus+1, parent) == p) ? antitheta : theta);
    }
    
    return tmp;
}

void SamplerRfunction::get_recombination_distribution(DescentGraph* dg, 
                                     unsigned person_id, 
                                     enum phased_trait parent_trait, 
                                     enum parentage parent, 
                                     double* dist) {
    
    switch(parent_trait) {
        case TRAIT_UU :
            dist[TRAIT_U] = 1.0;
            dist[TRAIT_A] = 0.0;
            return;
            
        case TRAIT_AU :
        case TRAIT_UA :
            break;
            
        case TRAIT_AA :
            dist[TRAIT_U] = 0.0;
            dist[TRAIT_A] = 1.0;
            return;
    }
    
    double tmp0 = 0.5;
    double tmp1 = 0.5;
    
    if(locus != 0) {
        bool cross = dg->get(person_id, locus-1, parent) != 0;
        
        tmp0 *= (cross ? theta2     : antitheta2);
        tmp1 *= (cross ? antitheta2 : theta2);
    }
    
    if(locus != (map->num_markers() - 1)) {
        bool cross = dg->get(person_id, locus+1, parent) != 0;
        
        tmp0 *= (cross ? theta     : antitheta);
        tmp1 *= (cross ? antitheta : theta);
    }
    
    double total = tmp0 + tmp1;
    
    if(parent_trait == TRAIT_AU) {
        dist[TRAIT_U] = tmp1 / total;
        dist[TRAIT_A] = 1.0 - dist[TRAIT_U];
    }
    else {
        dist[TRAIT_U] = tmp0 / total;
        dist[TRAIT_A] = 1.0 - dist[TRAIT_U];
    }
    
}

void SamplerRfunction::sample(vector<int>& pmk) {
    double prob_dist[NUM_ALLELES];
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        pmk[peel_id] = i;
        prob_dist[i] = pmatrix_presum.get(pmk);        
    }
    
    normalise(prob_dist);
        
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
        
        tmp *= transmission[0][transmission_index(mat_trait, pat_trait, kid_trait)];
        //tmp *= (get_recombination_probability(dg, peel_id, mat_trait, kid_trait, MATERNAL) *
        //        get_recombination_probability(dg, peel_id, pat_trait, kid_trait, PATERNAL));
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
    
    enum phased_trait kid_trait;
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
        
        for(unsigned int c = 0; c < children.size(); ++c) {
            unsigned int child_id = children[c];
            Person* child = ped->get_by_index(child_id);
            
            kid_trait = static_cast<enum phased_trait>(indices[pmatrix_index][child_id]);
            
            if(child->get_maternalid() == peel_id) {
                pat_trait = static_cast<enum phased_trait>(indices[pmatrix_index][child->get_paternalid()]);
            }
            else {
                mat_trait = static_cast<enum phased_trait>(indices[pmatrix_index][child->get_maternalid()]);
            }
            
            child_prob *= transmission[c][transmission_index(mat_trait, pat_trait, kid_trait)];
            //child_prob *= (get_recombination_probability(dg, child_id, mat_trait, kid_trait, MATERNAL) *
            //               get_recombination_probability(dg, child_id, pat_trait, kid_trait, PATERNAL));
            if(child_prob == 0.0)
                break;
        }
        
        tmp *= child_prob;
        
        pmatrix_presum.set(presum_index, tmp);
        
        total += tmp;
    }
    
    pmatrix.set(pmatrix_index, total);
}

void SamplerRfunction::setup_transmission_cache() {
        
    if(peel->get_type() == PARENT_PEEL) {
        for(unsigned int c = 0; c < peel->get_cutset_size(); ++c) {
            unsigned int child_id = peel->get_cutnode(c);
            Person* child = ped->get_by_index(child_id);
            
            if(not child->is_parent(peel_id))
                continue;
            
            transmission.push_back(new double[64]);
            children.push_back(child_id);
        }
    }
    else {
        transmission.push_back(new double[64]);
    }
}

void SamplerRfunction::teardown_transmission_cache() {
    
    for(unsigned int i = 0; i < transmission.size(); ++i) {
        delete[] transmission[i];
    }
    
    
    transmission.clear();
    children.clear();
}

void SamplerRfunction::preevaluate_init(DescentGraph* dg) {
    populate_transmission_cache(dg);
}

void SamplerRfunction::populate_transmission_cache(DescentGraph* dg) {

    if(peel->get_type() == PARENT_PEEL) {
        for(unsigned int i = 0; i < children.size(); ++i) {
            transmission_matrix(dg, children[i], transmission[i]);
        }
    }
    else {
        transmission_matrix(dg, peel_id, transmission[0]);
    }
}

unsigned int SamplerRfunction::transmission_index(enum phased_trait mat_trait, 
                                                  enum phased_trait pat_trait, 
                                                  enum phased_trait kid_trait) {

    return (mat_trait * 16) + (pat_trait * 4) + kid_trait;
}

void SamplerRfunction::transmission_matrix(DescentGraph* dg, int kid_id, double* tmatrix) {
    /*
    //#pragma omp parallel for
    for(int i = 0; i < 64; ++i) {
        int tmp = i;
        enum phased_trait mat_trait = static_cast<enum phased_trait>(tmp / 16); tmp %= 16;
        enum phased_trait pat_trait = static_cast<enum phased_trait>(tmp / 4);  tmp %= 4;
        enum phased_trait kid_trait = static_cast<enum phased_trait>(tmp);
    
        tmatrix[i] = get_recombination_probability(dg, kid_id, mat_trait, kid_trait, MATERNAL) *  \
                     get_recombination_probability(dg, kid_id, pat_trait, kid_trait, PATERNAL);
    }
    
    //#pragma omp parallel for
    for(int i = 0; i < 64; i += 4) {
        normalise(&tmatrix[i]);
    }
    */
    
    double mat_dist[2];
    double pat_dist[2];
    int mindex, pindex;
    
    for(int i = 0; i < 4; ++i) {
        mindex = 16 * i;
        get_recombination_distribution(dg, kid_id, static_cast<enum phased_trait>(i), MATERNAL, mat_dist);
        
        for(int j = 0; j < 4; ++j) {
            pindex = mindex + (4 * j);
            get_recombination_distribution(dg, kid_id, static_cast<enum phased_trait>(j), PATERNAL, pat_dist);
            
            // U = 0, A = 1
            // UU = 0, AU = 1, UA = 2, AA = 3
            tmatrix[pindex + TRAIT_UU] = mat_dist[TRAIT_U] * pat_dist[TRAIT_U];
            tmatrix[pindex + TRAIT_AU] = mat_dist[TRAIT_A] * pat_dist[TRAIT_U];
            tmatrix[pindex + TRAIT_UA] = mat_dist[TRAIT_U] * pat_dist[TRAIT_A];
            tmatrix[pindex + TRAIT_AA] = mat_dist[TRAIT_A] * pat_dist[TRAIT_A];
        }
    }
}

