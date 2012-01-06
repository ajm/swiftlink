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


Rfunction::Rfunction(PeelOperation po, Pedigree* p, GeneticMap* m, Rfunction* prev1, Rfunction* prev2) : 
      map(m),
      ped(p),
      offset(0.0),
      pmatrix(po.get_cutset_size(), NUM_ALLELES),
      pmatrix_presum(po.get_cutset_size() + 1, NUM_ALLELES),
      peel(po), 
      previous_rfunction1(prev1),
      previous_rfunction2(prev2),
      function_used(false),
      theta(0.0),
      antitheta(0.0),
      theta2(0.0),
      antitheta2(0.0),
      indices(peel.get_index_values()),
      size(pow(4.0, peel.get_cutset_size())) {

    pmatrix.set_keys(peel.get_cutset());
          
    // XXX temporary, neater way of doing this?
    vector<unsigned> tmp(peel.get_cutset());
    tmp.push_back(peel.get_peelnode());
    
    pmatrix_presum.set_keys(tmp);
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
    function_used(r.function_used),
    theta(r.theta),
    antitheta(r.antitheta),
    theta2(r.theta2),
    antitheta2(r.antitheta2),
    indices(r.indices),
    size(r.size) {}
    
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
        function_used = rhs.function_used;
        theta = rhs.theta;
        antitheta = rhs.antitheta;
        theta2 = rhs.theta2;
        antitheta2 = rhs.antitheta2;
        indices = rhs.indices;
    }
    
    return *this;
}

// this function is the same for Traits and Sampling
void Rfunction::evaluate_partner_peel(unsigned int pmatrix_index, unsigned locus) {
    double tmp = 0.0;
    double total = 0.0;
    
    unsigned int offset = 1 << (2 * peel.get_cutset_size());
    unsigned int presum_index;
    
    enum phased_trait partner_trait;
    unsigned partner_id = peel.get_peelnode();
        
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        partner_trait = static_cast<enum phased_trait>(i);        
        presum_index = pmatrix_index + (offset * i);
        
        (*indices)[pmatrix_index][partner_id] = i;
        
        tmp = get_trait_probability(partner_id, partner_trait, locus);
        
        if(tmp == 0.0)
            continue;
        
        tmp *= (previous_rfunction1 == NULL ? 1.0 : previous_rfunction1->get((*indices)[pmatrix_index])) * \
               (previous_rfunction2 == NULL ? 1.0 : previous_rfunction2->get((*indices)[pmatrix_index]));
        
        pmatrix_presum.set(presum_index, tmp);
        
        total += tmp;
    }
    
    pmatrix.set(pmatrix_index, total);
}

void Rfunction::evaluate_element(
                    unsigned int pmatrix_index, 
                    DescentGraph* dg, 
                    unsigned locus) {
    
    // XXX could remove this with some inheritance?
    // RfunctionChild RfunctionParent?
    switch(peel.get_type()) {
        
        case CHILD_PEEL :
            evaluate_child_peel(pmatrix_index, dg, locus);
            break;
            
        case PARTNER_PEEL :
        case LAST_PEEL :
            evaluate_partner_peel(pmatrix_index, locus);
            break;
        
        case PARENT_PEEL :
            evaluate_parent_peel(pmatrix_index, dg, locus);
            break;
        
        default :
            fprintf(stderr, "error: default should never be reached! (%s:%d)\n", __FILE__, __LINE__);
            abort();
    }
}

bool Rfunction::legal_genotype(unsigned personid, unsigned locus, enum phased_trait g) {

    //return true;

    Person* p = ped->get_by_index(personid);
    
    return p->legal_genotype(locus, g);
}

void Rfunction::evaluate(DescentGraph* dg, unsigned locus, double offset) {
    //PeelMatrixKey k(ped->num_members());
    //vector<unsigned> q;
    //unsigned ndim = peel.get_cutset_size();
    //unsigned tmp;
    //bool skip = false;
    
    pmatrix.reset();
    pmatrix_presum.reset();
    
    // crucial for TraitRfunction
    this->offset = offset;
    
    if(locus != (map->num_markers() - 1)) {
        theta = map->get_theta(locus);
        antitheta = map->get_inversetheta(locus);
    }
    
    if(locus != 0) {
        theta2 = map->get_theta(locus-1);
        antitheta2 = map->get_inversetheta(locus-1);
    }
    
    /*
    // nothing in the cutset to be enumerated
    if(peel.get_type() == LAST_PEEL) {
        evaluate_partner_peel(k, locus);
        return;
    }
    */
    
    for(unsigned int i = 0; i < size; ++i) {
        //fprintf(stderr, "%d\n", i);
        evaluate_element(i, dg, locus);
    }
    //fprintf(stderr, "====\n");
}

