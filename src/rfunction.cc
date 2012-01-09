using namespace std;

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


Rfunction::Rfunction(Pedigree* p, GeneticMap* m, unsigned int locus, PeelOperation* po, Rfunction* prev1, Rfunction* prev2) : 
    map(m),
    ped(p),
    offset(0.0),
    pmatrix(po->get_cutset_size(), NUM_ALLELES),
    pmatrix_presum(po->get_cutset_size() + 1, NUM_ALLELES),
    peel(po), 
    previous_rfunction1(prev1),
    previous_rfunction2(prev2),
    locus(locus),
    theta(0.0),
    antitheta(0.0),
    theta2(0.0),
    antitheta2(0.0),
    indices(peel->get_index_values()),
    index_offset(1 << (2 * peel->get_cutset_size())),
    size(pow((double)NUM_ALLELES, (int)peel->get_cutset_size())),
    peel_id(peel->get_peelnode()) {
    
    pmatrix.set_keys(peel->get_cutset());
          
    // XXX temporary, neater way of doing this?
    vector<unsigned> tmp(peel->get_cutset());
    tmp.push_back(peel->get_peelnode());
    
    pmatrix_presum.set_keys(tmp);
    
    
    if(locus != (map->num_markers() - 1)) {
        theta = map->get_theta(locus);
        antitheta = map->get_inversetheta(locus);
    }
    
    if(locus != 0) {
        theta2 = map->get_theta(locus-1);
        antitheta2 = map->get_inversetheta(locus-1);
    }        
}

Rfunction::Rfunction(const Rfunction& r) :
    map(r.map),
    ped(r.ped),
    offset(r.offset),
    pmatrix(r.pmatrix),
    pmatrix_presum(r.pmatrix_presum),
    peel(r.peel),
    previous_rfunction1(r.previous_rfunction1),
    previous_rfunction2(r.previous_rfunction2),
    locus(r.locus),
    theta(r.theta),
    antitheta(r.antitheta),
    theta2(r.theta2),
    antitheta2(r.antitheta2),
    indices(r.indices),
    index_offset(r.index_offset),
    size(r.size),
    peel_id(r.peel_id) {}
    
Rfunction& Rfunction::operator=(const Rfunction& rhs) {

    if(&rhs != this) {
        pmatrix = rhs.pmatrix;
        pmatrix_presum = rhs.pmatrix_presum;
        peel = rhs.peel;
        map = rhs.map;
        ped = rhs.ped;
        offset = rhs.offset;
        previous_rfunction1 = rhs.previous_rfunction1;
        previous_rfunction2 = rhs.previous_rfunction2;
        locus = rhs.locus;
        theta = rhs.theta;
        antitheta = rhs.antitheta;
        theta2 = rhs.theta2;
        antitheta2 = rhs.antitheta2;
        indices = rhs.indices;
        index_offset = rhs.index_offset;
        size = rhs.size;
        peel_id = rhs.peel_id;
    }
    
    return *this;
}

// this function is the same for Traits and Sampling
void Rfunction::evaluate_partner_peel(unsigned int pmatrix_index) {
    double tmp = 0.0;
    double total = 0.0;
    
    unsigned int presum_index;
    
    enum phased_trait partner_trait;
    
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        partner_trait = static_cast<enum phased_trait>(i);        
        presum_index = pmatrix_index + (index_offset * i);
        
        (*indices)[pmatrix_index][peel_id] = i;
        
        tmp = get_trait_probability(peel_id, partner_trait);
        
        if(tmp == 0.0)
            continue;
        
        tmp *= (previous_rfunction1 == NULL ? 1.0 : previous_rfunction1->get((*indices)[pmatrix_index])) * \
               (previous_rfunction2 == NULL ? 1.0 : previous_rfunction2->get((*indices)[pmatrix_index]));
        
        pmatrix_presum.set(presum_index, tmp);
        
        total += tmp;
    }
    
    pmatrix.set(pmatrix_index, total);
}

void Rfunction::evaluate_element(unsigned int pmatrix_index, DescentGraph* dg) {
    // XXX could remove this with some inheritance?
    // RfunctionChild RfunctionParent?
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

void Rfunction::evaluate(DescentGraph* dg, double offset) {
    pmatrix.reset();
    pmatrix_presum.reset();
    
    // crucial for TraitRfunction
    this->offset = offset;
    
    for(unsigned int i = 0; i < size; ++i) {
        evaluate_element(i, dg);
    }
}
