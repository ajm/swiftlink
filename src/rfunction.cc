#include <cmath>
#include <vector>

#include "rfunction.h"
#include "peeling.h"
#include "peel_matrix.h"
#include "genotype.h"
#include "pedigree.h"
#include "descent_graph.h"
#include "trait.h"
#include "genetic_map.h"

using namespace std;


Rfunction::Rfunction(Pedigree* p, GeneticMap* m, unsigned int locus, PeelOperation* po, vector<Rfunction*> previous, bool sex_linked) :
    map(m),
    ped(p),
    offset(0),
    pmatrix(po->get_cutset_size(), NUM_ALLELES),
    pmatrix_presum(po->get_cutset_size() + 1, NUM_ALLELES),
    peel(po), 
    previous_rfunctions(previous),
    locus(locus),
    indices(peel->get_index_values()),
    valid_indices(peel->get_matrix_indices(locus)),
    valid_lod_indices(peel->get_lod_indices()),
    index_offset(1 << (2 * peel->get_cutset_size())),
    size(pow((double)NUM_ALLELES, (int)peel->get_cutset_size())),
    peel_id(peel->get_peelnode()),
    theta(0.0),
    antitheta(1.0),
    theta2(0.0),
    antitheta2(1.0),
    sex_linked(sex_linked) {
    
    pmatrix.set_keys(peel->get_cutset());
          
    // XXX temporary, neater way of doing this?
    vector<unsigned> tmp(peel->get_cutset());
    tmp.push_back(peel->get_peelnode());
    
    pmatrix_presum.set_keys(tmp);      
}

Rfunction::Rfunction(const Rfunction& rhs) :
    map(rhs.map),
    ped(rhs.ped),
    offset(rhs.offset),
    pmatrix(rhs.pmatrix),
    pmatrix_presum(rhs.pmatrix_presum),
    peel(rhs.peel),
    previous_rfunctions(rhs.previous_rfunctions),
    locus(rhs.locus),
    indices(rhs.indices),
    valid_indices(rhs.valid_indices),
    valid_lod_indices(rhs.valid_lod_indices),
    index_offset(rhs.index_offset),
    size(rhs.size),
    peel_id(rhs.peel_id),
    theta(rhs.theta),
    antitheta(rhs.antitheta),
    theta2(rhs.theta2),
    antitheta2(rhs.antitheta2),
    sex_linked(rhs.sex_linked) {}
    
Rfunction& Rfunction::operator=(const Rfunction& rhs) {

    if(&rhs != this) {
        pmatrix = rhs.pmatrix;
        pmatrix_presum = rhs.pmatrix_presum;
        peel = rhs.peel;
        map = rhs.map;
        ped = rhs.ped;
        offset = rhs.offset;
        previous_rfunctions = rhs.previous_rfunctions;
        locus = rhs.locus;
        indices = rhs.indices;
        valid_indices = rhs.valid_indices;
        valid_lod_indices = rhs.valid_lod_indices;
        index_offset = rhs.index_offset;
        size = rhs.size;
        peel_id = rhs.peel_id;
        theta = rhs.theta;
        antitheta = rhs.antitheta;
        theta2 = rhs.theta2;
        antitheta2 = rhs.antitheta2;
        sex_linked = rhs.sex_linked;
    }
    
    return *this;
}

enum phased_trait Rfunction::get_phased_trait(enum phased_trait m, enum phased_trait p, 
                                                   int maternal_allele, int paternal_allele, enum sex child_sex) {
                                                   
    bool m_affected = affected_trait(m, maternal_allele);
    bool p_affected = affected_trait(p, paternal_allele);
    enum phased_trait pt;
    
    if(sex_linked and child_sex == MALE) {
        return m_affected ? TRAIT_AA : TRAIT_UU;
    }

    if(m_affected) {
        pt = p_affected ? TRAIT_AA : TRAIT_AU;
    }
    else {
        pt = p_affected ? TRAIT_UA : TRAIT_UU;
    }
    
    return pt;
}

// this function is the same for Traits and Sampling
void Rfunction::evaluate_partner_peel(unsigned int pmatrix_index) {
    double tmp = 0.0;
    double total = 0.0;
    
    unsigned int presum_index;
    
    enum phased_trait partner_trait;
    
    
    for(unsigned i = 0; i < 4; ++i) {
        partner_trait = static_cast<enum phased_trait>(i);        
        presum_index = pmatrix_index + (index_offset * i);
        
        indices[pmatrix_index][peel_id] = i;
        
        tmp = trait_cache[i];
        if(tmp == 0.0)
            continue;
        
        for(unsigned int j = 0; j < previous_rfunctions.size(); ++j) {
            tmp *= previous_rfunctions[j]->get(indices[pmatrix_index]);
        }
        
        pmatrix_presum.set(presum_index, tmp);
        
        total += tmp;
    }
    
    pmatrix.set(pmatrix_index, total);
/*
    if(total == 0.0) {
        fprintf(stderr, "ERROR: evaluate partner peel = ZERO\n");
    }
*/
}

void Rfunction::evaluate_element(unsigned int pmatrix_index, DescentGraph* dg) {
    switch(peel->get_type()) {
        
        case CHILD_PEEL :
            evaluate_child_peel(pmatrix_index, dg);
            break;
            
        case PARTNER_PEEL :
        case LAST_PEEL :
            evaluate_partner_peel(pmatrix_index);
            break;
        
        case PARENT_PEEL :
            evaluate_parent_peel(pmatrix_index, dg);
            break;
        
        default :
            fprintf(stderr, "error: default should never be reached! (%s:%d)\n", __FILE__, __LINE__);
            abort();
    }
}

bool Rfunction::legal_genotype(unsigned personid, enum phased_trait g) {
    Person* p = ped->get_by_index(personid);
    
    return p->legal_genotype(locus, g);
}

void Rfunction::evaluate(DescentGraph* dg, unsigned int offset) {
    //pmatrix.reset();
    //pmatrix_presum.reset();
    
    // transmission probability cache
    preevaluate_init(dg);
    
    // crucial for TraitRfunction
    this->offset = offset;
    
    //#pragma omp parallel for
    
    // calculate lod score
    if(offset != 0) {
        //for(unsigned int i = 0; i < size; ++i) {
        for(unsigned int i = 0; i < valid_lod_indices->size(); ++i) {
            evaluate_element((*valid_lod_indices)[i], dg);
        }
    }
    // running locus sampler
    else {
        // this is only for the SamplerRfunction at the moment
        for(unsigned int i = 0; i < valid_indices->size(); ++i) {
            evaluate_element((*valid_indices)[i], dg);
        }
    }
}

void Rfunction::normalise(double* p) {
    double total = p[0] + p[1] + p[2] + p[3];
    
    if(total == 0.0)
        return;
    
    for(int i = 0; i < 4; ++i) {
        p[i] /= total;
    }
}

