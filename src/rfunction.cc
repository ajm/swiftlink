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
      antitheta2(0.0) {

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
    antitheta2(r.antitheta2) {}
    
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
    }
    
    return *this;
}

bool Rfunction::contains_node(unsigned node) {
    vector<unsigned>& cutset = peel.get_cutset();
    
    return find(cutset.begin(), cutset.end(), node) != cutset.end();
}

bool Rfunction::contains_cutnodes(vector<unsigned>& nodes) {
    vector<unsigned>& cutset = peel.get_cutset();
    
    for(unsigned i = 0; i < nodes.size(); ++i) {
        if(find(cutset.begin(), cutset.end(), nodes[i]) == cutset.end())
            return false;
    }
    
    return true;
}

bool Rfunction::is_used() {
    return function_used;
}

void Rfunction::set_used() {
    function_used = true;
}

void Rfunction::generate_key(PeelMatrixKey& pmatrix_index, vector<unsigned int>& assignments) {
    pmatrix_index.reassign(peel.get_cutset(), assignments);
}

void Rfunction::summation(PeelMatrixKey& pmatrix_index, unsigned person_id) {
    double tmp = 0.0;
    enum phased_trait pt;
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        pt = static_cast<enum phased_trait>(i);
        pmatrix_index.add(person_id, pt);
        
        tmp += pmatrix_presum.get(pmatrix_index);
    }
    
    pmatrix.set(pmatrix_index, tmp);
}

/*
enum trait Rfunction::get_trait(enum phased_trait p, enum parentage parent) {
    
    switch(parent) {
        case MATERNAL:
            return (((p == TRAIT_UU) or (p == TRAIT_UA)) ? TRAIT_U : TRAIT_A);
        
        case PATERNAL:
            return (((p == TRAIT_UU) or (p == TRAIT_AU)) ? TRAIT_U : TRAIT_A);
        
        default:
            break;
    }
    
    abort();
}
*/

void Rfunction::get_traits(enum phased_trait p, enum trait& mat, enum trait& pat) {
    switch(p) {
        case TRAIT_UU:
            mat = pat = TRAIT_U;
            break;
        case TRAIT_AU:
            mat = TRAIT_A;
            pat = TRAIT_U;
            break;
        case TRAIT_UA:
            mat = TRAIT_U;
            pat = TRAIT_A;
            break;
        case TRAIT_AA:
            mat = pat = TRAIT_A;
            break;
    }
}

// this function is the same for Traits and Sampling
void Rfunction::evaluate_partner_peel(PeelMatrixKey& pmatrix_index, unsigned locus) {
    
    double tmp = 0.0;
    double total = 0.0;
    unsigned partner_id;
    enum phased_trait partner_trait;    
    
    partner_id = peel.get_peelnode();
        
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        partner_trait = static_cast<enum phased_trait>(i);
        pmatrix_index.add(partner_id, partner_trait);
        
        tmp = get_trait_probability(partner_id, partner_trait, locus) * \
            (previous_rfunction1 == NULL ? 1.0 : previous_rfunction1->get(pmatrix_index)) * \
            (previous_rfunction2 == NULL ? 1.0 : previous_rfunction2->get(pmatrix_index));
        
        pmatrix_presum.set(pmatrix_index, tmp);
        total += tmp;
    }
    
    //summation(pmatrix_index, partner_id);
    pmatrix.set(pmatrix_index, total);
}

void Rfunction::evaluate_element(
                    PeelMatrixKey& pmatrix_index, 
                    DescentGraph* dg, 
                    unsigned locus) {
    
    // XXX could remove this with some inheritance?
    // RfunctionChild RfunctionParent?
    switch(peel.get_type()) {
        
        case CHILD_PEEL :
            evaluate_child_peel(pmatrix_index, dg, locus);
            break;
            
        case PARTNER_PEEL :
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

void Rfunction::evaluate(DescentGraph* dg, unsigned locus, double offset) {
    PeelMatrixKey k(ped->num_members());
    vector<unsigned> q;
    unsigned ndim = peel.get_cutset_size();
    unsigned tmp;
    
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
    
    
    // nothing in the cutset to be enumerated
    if(peel.get_type() == LAST_PEEL) {
        evaluate_partner_peel(k, locus);
        return;
    }
    
    // generate all assignments to iterate through n-dimensional matrix
    
    // initialise to the first element of matrix
    for(unsigned i = 0; i < ndim; ++i) {
        q.push_back(0);
    }

    // enumerate all elements in ndim-dimenstional matrix
    while(not q.empty()) {
        
        if(q.size() == ndim) {
            generate_key(k, q);
            evaluate_element(k, dg, locus);
        }
        
        tmp = q.back() + 1;
        q.pop_back();
        
        if(tmp < NUM_ALLELES) {
            q.push_back(tmp);
            tmp = ndim - q.size();
            // fill out rest with zeroes
            for(unsigned i = 0; i < tmp; ++i) {
                q.push_back(0);
            }
        }
    }
}
