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
      temperature(0.0),
      pmatrix(po.get_cutset_size(), NUM_ALLELES),
      pmatrix_presum(po.get_cutset_size() + 1, NUM_ALLELES),
      peel(po), 
      previous_rfunction1(prev1),
      previous_rfunction2(prev2),
      function_used(false) {

    pmatrix.set_keys(peel.get_cutset());
          
    // very temporary, just messing around... XXX
    vector<unsigned> tmp(peel.get_cutset());
    tmp.push_back(peel.get_peelnode());
    
    pmatrix_presum.set_keys(tmp);
    
          /*
    printf("prev1 = %s\nprev2 = %s\n\n", \
            previous_rfunction1 == NULL ? "NULL" : "NOT NULL", \
            previous_rfunction2 == NULL ? "NULL" : "NOT NULL");
          */
}

Rfunction::Rfunction(const Rfunction& r) :
    map(r.map),
    ped(r.ped),
    offset(r.offset),
    temperature(r.temperature),
    pmatrix(r.pmatrix),
    pmatrix_presum(r.pmatrix_presum),
    peel(r.peel),
    previous_rfunction1(r.previous_rfunction1),
    previous_rfunction2(r.previous_rfunction2),
    function_used(r.function_used) {}
    
Rfunction& Rfunction::operator=(const Rfunction& rhs) {

    if(&rhs != this) {
        pmatrix = rhs.pmatrix;
        pmatrix_presum = rhs.pmatrix_presum;
        peel = rhs.peel;
        map = rhs.map;
        ped = rhs.ped;
        offset = rhs.offset;
        temperature = rhs.temperature;
        previous_rfunction1 = rhs.previous_rfunction1;
        previous_rfunction2 = rhs.previous_rfunction2;
        function_used = rhs.function_used;
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

bool Rfunction::affected_trait(enum phased_trait pt, int allele) {
    
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

enum phased_trait Rfunction::get_phased_trait(
                    enum phased_trait m, enum phased_trait p, 
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

void Rfunction::evaluate_child_peel(
                    PeelMatrixKey& pmatrix_index, 
                    DescentGraph* dg,
                    unsigned locus) {
    
    double recombination_prob;
    double transmission_prob;
    double disease_prob;
    double old_prob1;
    double old_prob2;
    double tmp = 0.0;
    
    enum phased_trait kid_trait;
    enum phased_trait mat_trait;
    enum phased_trait pat_trait;
    
    unsigned kid_id = peel.get_peelnode();
    Person* kid = ped->get_by_index(kid_id);
    
    mat_trait = pmatrix_index.get(kid->get_maternalid());
    pat_trait = pmatrix_index.get(kid->get_paternalid());
    
    
    for(unsigned i = 0; i < NUM_ALLELES; ++i) {
        kid_trait = static_cast<enum phased_trait>(i);
        pmatrix_index.add(kid_id, kid_trait);
        pmatrix_presum.set(pmatrix_index, 0.0);
    }
    
    
    // iterate over all descent graphs to determine child trait 
    // based on parents' traits
    for(int i = 0; i < 2; ++i) {        // maternal
        for(int j = 0; j < 2; ++j) {    // paternal
            
            kid_trait = get_phased_trait(mat_trait, pat_trait, i, j);
            
            pmatrix_index.add(kid_id, kid_trait);
            
            disease_prob        = get_trait_probability(kid_id, kid_trait, locus);
            transmission_prob   = get_transmission_probability(mat_trait) *  \
                                  get_transmission_probability(pat_trait);
            recombination_prob  = !dg ? 1.0 : get_recombination_probability(dg, locus, kid_id, i, j);
            old_prob1           = previous_rfunction1 != NULL ? previous_rfunction1->get(pmatrix_index) : 1.0;
            old_prob2           = previous_rfunction2 != NULL ? previous_rfunction2->get(pmatrix_index) : 1.0;
            
            tmp = (disease_prob * transmission_prob * recombination_prob * old_prob1 * old_prob2);
            
            pmatrix_presum.add(pmatrix_index, tmp);
        }
    }
    
    summation(pmatrix_index, kid_id);
}

// TODO XXX this is a mess
//
void Rfunction::evaluate_parent_peel(
                    PeelMatrixKey& pmatrix_index, 
                    DescentGraph* dg,
                    unsigned locus) {
    
    double disease_prob;
    double recombination_prob;
    double old_prob1;
    double old_prob2;
    
    
    unsigned parent_id = peel.get_peelnode();
    unsigned child_node = peel.get_cutnode(0);
    
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
    unsigned other_parent_id = ismother ? \
    p->get_paternalid() : \
    p->get_maternalid();
    enum phased_trait pivot_trait;
    enum phased_trait parent_trait;
    enum phased_trait other_trait;
    double tmp = 0.0;
    
    other_trait = pmatrix_index.get(other_parent_id);
    
    
    for(unsigned a = 0; a < NUM_ALLELES; ++a) {
        parent_trait = static_cast<enum phased_trait>(a);
        pmatrix_index.add(parent_id, parent_trait);
        
        disease_prob = get_trait_probability(parent_id, parent_trait, locus);
        
        old_prob1 = previous_rfunction1 != NULL ? previous_rfunction1->get(pmatrix_index) : 1.0;
        old_prob2 = previous_rfunction2 != NULL ? previous_rfunction2->get(pmatrix_index) : 1.0;
        
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
                    
                    if(pivot_trait != pmatrix_index.get(child->get_internalid()))
                        continue;
                    
                    recombination_prob = !dg ? 0.25 : get_recombination_probability(dg, locus, child->get_internalid(), i, j);
                    
                    child_tmp += recombination_prob; //(disease_prob * recombination_prob * old_prob1 * old_prob2);
                }
            }
            
            child_prob *= child_tmp;
        }
        
        //tmp += (child_prob * disease_prob * old_prob1 * old_prob2);
        tmp = (child_prob * disease_prob * old_prob1 * old_prob2);
    
        pmatrix_presum.set(pmatrix_index, tmp);
    }
    
    summation(pmatrix_index, parent_id);
}

void Rfunction::evaluate_partner_peel(PeelMatrixKey& pmatrix_index, unsigned locus) {
    
    double tmp = 0.0;
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
    }
    
    summation(pmatrix_index, partner_id);
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

void Rfunction::evaluate(DescentGraph* dg, unsigned locus, double offset, double temperature) {
    PeelMatrixKey k;
    vector<unsigned> q;
    unsigned ndim = peel.get_cutset_size();
    unsigned tmp;
    
    // crucial for TraitRfunction
    this->offset = offset;
    
    // necessary for SamplerRfunction
    this->temperature = temperature;
    
    
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

