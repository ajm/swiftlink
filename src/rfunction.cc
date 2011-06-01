using namespace std;

#include <cmath>

#include "rfunction.h"
#include "peeling.h"
#include "peel_matrix.h"
#include "genotype.h"
#include "pedigree.h"
#include "descent_graph.h"
#include "trait.h"
#include "genetic_map.h"


Rfunction::Rfunction(PeelOperation po, Pedigree* p, GeneticMap* m, unsigned alleles, vector<Rfunction>& previous_functions, unsigned index)
    : pmatrix(po.get_cutset_size(), alleles), 
      peel(po), 
      num_alleles(alleles), 
      map(m),
      ped(p),
      function_used(false),
      function_index(index) {
        
    if(alleles != 4)
        abort(); // XXX assumption for now...

    pmatrix.set_keys(peel.get_cutset());
    
    find_previous_functions(previous_functions);

    printf("RFUNCTION: %d ", function_index);
    peel.print();
    printf("(deps: %d %d)\n", 
        previous_rfunction1 == NULL ? -1 : previous_rfunction1->function_index, 
        previous_rfunction2 == NULL ? -1 : previous_rfunction2->function_index);
}

void Rfunction::find_previous_functions(vector<Rfunction>& functions) {
    switch(peel.get_type()) {
        case CHILD_PEEL:
            find_child_functions(functions);
            break;
            
        case PARTNER_PEEL:
            find_partner_functions(functions);
            break;
            
        case PARENT_PEEL:
            find_parent_functions(functions);
            break;
            
        case LAST_PEEL:
            find_last_functions(functions);
            break;
            
        default:
            fprintf(stderr, "error: default should never be reached! (%s:%d)\n", __FILE__, __LINE__);
            abort();
    }
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

void Rfunction::find_function_containing(vector<Rfunction>& functions, vector<unsigned>& nodes, Rfunction** func) {
    
    for(unsigned i = 0; i < functions.size(); ++i) {
        if(functions[i].is_used()) {
            continue;
        }
        
        if(functions[i].contains_cutnodes(nodes)) {
            *func = &functions[i];
            functions[i].set_used();
            break;
        }
    }
}

void Rfunction::find_child_functions(vector<Rfunction>& functions) {

    vector<unsigned>& peelset = peel.get_peelset();
    
    Person* p = ped->get_by_index(peelset[0]);

    vector<unsigned> tmp;
    tmp.push_back(p->get_maternalid());
    tmp.push_back(p->get_paternalid());
    
    previous_rfunction1 = NULL;
    previous_rfunction2 = NULL;
    
    find_function_containing(functions, tmp, &previous_rfunction1);
    
    // don't even bother looking if the child is a leaf
    if(p->isleaf()) {
        return;
    }
    
    find_function_containing(functions, peelset, &previous_rfunction2);
    
    if(previous_rfunction2 == NULL) {
        fprintf(stderr, "error: non-leaf child node not found in any previous function (%s:%d)\n", 
            __FILE__, __LINE__);
        abort();
    }
}

void Rfunction::find_partner_functions(vector<Rfunction>& functions) {
    
    vector<unsigned>& peelset = peel.get_peelset();
    vector<unsigned>& cutset = peel.get_cutset();
    
    Person* p = ped->get_by_index(peelset[0]);
    
    vector<unsigned> tmp;
    tmp.push_back(peelset[0]);
    
    for(unsigned i = 0; i < p->num_mates(); ++i) {
        unsigned tmp2 = p->get_mate(i)->get_internalid();
        for(unsigned j = 0; j < cutset.size(); ++j) {
            if(cutset[j] == tmp2) {
                tmp.push_back(tmp2);
            }
        }
    }
        
    if(tmp.size() != 2) {
        fprintf(stderr, "error: ambiguous which mate is being peeled on to (%s:%d)\n", 
            __FILE__, __LINE__);
        abort();
    }
    
    previous_rfunction1 = NULL;
    previous_rfunction2 = NULL;
    
    find_function_containing(functions, tmp, &previous_rfunction1);
    
    
    if(p->isfounder()) {
        return;
    }
    
    find_function_containing(functions, peelset, &previous_rfunction2);
    
    if(previous_rfunction2 == NULL) {
        fprintf(stderr, "error: non-founder node not found in any previous function (%s:%d)\n", 
            __FILE__, __LINE__);
        abort();
    }
}

void Rfunction::find_parent_functions(vector<Rfunction>& functions) {

    vector<unsigned>& peelset = peel.get_peelset();
    Person* mother;
    Person* father;
    Person* tmp;
    
    previous_rfunction1 = NULL;
    previous_rfunction2 = NULL;
    
    tmp = ped->get_by_index(peelset[0]);
    if(tmp->ismale()) {
        father = tmp;
        mother = ped->get_by_index(peelset[1]);
    }
    else {
        mother = tmp;
        father = ped->get_by_index(peelset[1]);
    }
    
    vector<unsigned> tmp2;
    tmp2.push_back(mother->get_internalid());
    
    find_function_containing(functions, tmp2, &previous_rfunction1);
    
    tmp2.clear();
    tmp2.push_back(father->get_internalid());
    
    find_function_containing(functions, tmp2, &previous_rfunction2);

    if((not mother->isfounder()) and (previous_rfunction1 == NULL)) {
        fprintf(stderr, "error: non-founder node not found in any previous function (%s:%d)\n", 
            __FILE__, __LINE__);
        abort();
    }

    if((not father->isfounder()) and (previous_rfunction2 == NULL)) {
        fprintf(stderr, "error: non-founder node not found in any previous function (%s:%d)\n", 
            __FILE__, __LINE__);
        abort();
    }
}

void Rfunction::find_last_functions(vector<Rfunction>& functions) {
    vector<unsigned>& peelset = peel.get_peelset();
    
    previous_rfunction1 = NULL;
    previous_rfunction2 = NULL;

    find_function_containing(functions, peelset, &previous_rfunction1);
    
    if(previous_rfunction1 == NULL) {
        fprintf(stderr, "error: last node not found in any previous function (%s:%d)\n", 
            __FILE__, __LINE__);
        abort();
    }
    
    find_function_containing(functions, peelset, &previous_rfunction2);
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

double Rfunction::get_disease_probability(unsigned person_id, enum phased_trait pt) {
    Person* per;
    
    per = ped->get_by_index(person_id);
    
    return per->get_disease_prob(pt);
}

// XXX this needs to be somewhere that is not at the same position as 
// a genetic marker
// TODO XXX make this generic so it can do arbitrary points between markers
double Rfunction::get_recombination_probability(
                    DescentGraph* dg, unsigned int locus_index, unsigned person_id,
                    int maternal_allele, int paternal_allele) {

    double tmp = 1.0;
    double half_recomb_prob;
    
    half_recomb_prob = map->get_theta_halfway(locus_index);
    
    tmp *= dg->get(person_id, locus_index,   MATERNAL) == maternal_allele ? 1.0 - half_recomb_prob : half_recomb_prob;
    tmp *= dg->get(person_id, locus_index+1, MATERNAL) == maternal_allele ? 1.0 - half_recomb_prob : half_recomb_prob;
            
    tmp *= dg->get(person_id, locus_index,   PATERNAL) == paternal_allele ? 1.0 - half_recomb_prob : half_recomb_prob;
    tmp *= dg->get(person_id, locus_index+1, PATERNAL) == paternal_allele ? 1.0 - half_recomb_prob : half_recomb_prob;
    
    return tmp;
}

void Rfunction::evaluate_last_peel(PeelMatrixKey& pmatrix_index) {  
    
    double tmp = 0.0;
    enum phased_trait last_trait;
    unsigned last_id;
    
    last_id = peel.get_peelnode(0);
            
    for(unsigned i = 0; i < num_alleles; ++i) {
        last_trait = static_cast<enum phased_trait>(i);
        
        pmatrix_index.add(last_id, last_trait);
        
        tmp += (\
                get_disease_probability(last_id, last_trait) * \
                previous_rfunction1->get(pmatrix_index) * \
                (previous_rfunction2 == NULL ? 1.0 : previous_rfunction2->get(pmatrix_index)) \
            );
    }
    
    pmatrix_index.add(last_id, (enum phased_trait) 0);
    pmatrix.set(pmatrix_index, tmp);
}

void Rfunction::evaluate_child_peel(
                    PeelMatrixKey& pmatrix_index, 
                    DescentGraph* dg,
                    unsigned locus_index) {
    
    double recombination_prob;
    double disease_prob;
    double old_prob1;
    double old_prob2;
    double tmp = 0.0;
    
    enum phased_trait kid_trait;
    enum phased_trait mat_trait;
    enum phased_trait pat_trait;
    
    unsigned kid_id = peel.get_peelnode(0);
    Person* kid = ped->get_by_index(kid_id);
    
    mat_trait = pmatrix_index.get(kid->get_maternalid());
    pat_trait = pmatrix_index.get(kid->get_paternalid());
    
    // iterate over all descent graphs to determine child trait 
    // based on parents' traits
    for(int i = 0; i < 2; ++i) {        // maternal
        for(int j = 0; j < 2; ++j) {    // paternal
            
            kid_trait = get_phased_trait(mat_trait, pat_trait, i, j);
            
            pmatrix_index.add(kid_id, kid_trait);
            
            disease_prob        = get_disease_probability(kid_id, kid_trait);
            recombination_prob  = !dg ? 0.25 : 0.25 * get_recombination_probability(dg, locus_index, kid_id, i, j);
            old_prob1           = previous_rfunction1 != NULL ? previous_rfunction1->get(pmatrix_index) : 1.0;
            old_prob2           = previous_rfunction2 != NULL ? previous_rfunction2->get(pmatrix_index) : 1.0;
            
            tmp += (disease_prob * recombination_prob * old_prob1 * old_prob2);
        }
    }
    
    pmatrix.set(pmatrix_index, tmp);
}

void Rfunction::evaluate_parent_peel(
                    PeelMatrixKey& pmatrix_index, 
                    DescentGraph* dg,
                    unsigned int locus_index) {
    
    double child_prob[4];
    double maternal_disease_prob;
    double paternal_disease_prob;
    double recombination_prob;
    double old_prob1;
    double old_prob2;
    
    enum phased_trait pivot_trait;
    unsigned piv_id = peel.get_cutnode(0);
    Person* p = ped->get_by_index(piv_id);
    unsigned mat_id = p->get_maternalid();
    unsigned pat_id = p->get_paternalid();
    
        
    for(unsigned i = 0; i < num_alleles; ++i)
        child_prob[i] = 0.0;
    
    for(unsigned mi = 0; mi < num_alleles; ++mi) {        // mother's genotype
        for(unsigned pi = 0; pi < num_alleles; ++pi) {    // father's genotype
            enum phased_trait m = static_cast<enum phased_trait>(mi);
            enum phased_trait p = static_cast<enum phased_trait>(pi);
            
            for(int i = 0; i < 2; ++i) {        // maternal inheritance
                for(int j = 0; j < 2; ++j) {    // paternal inheritance
                    pivot_trait = get_phased_trait(m, p, i, j);
                    
                    pmatrix_index.add(mat_id, m); // add mother
                    pmatrix_index.add(pat_id, p); // add father

                    maternal_disease_prob = get_disease_probability(mat_id, m);
                    paternal_disease_prob = get_disease_probability(pat_id, p);
                    recombination_prob    = !dg ? 0.25 : 0.25 * get_recombination_probability(dg, locus_index, piv_id, i, j);
                    old_prob1             = previous_rfunction1 != NULL ? previous_rfunction1->get(pmatrix_index) : 1.0;
                    old_prob2             = previous_rfunction2 != NULL ? previous_rfunction2->get(pmatrix_index) : 1.0;
                    
                    child_prob[pivot_trait] += \
                        (maternal_disease_prob * \
                         paternal_disease_prob * \
                         recombination_prob * \
                         old_prob1 * \
                         old_prob2 );
                }
            }
        }
    }
    
    for(unsigned i = 0; i < num_alleles; ++i) {
        pivot_trait = static_cast<enum phased_trait>(i);
        pmatrix_index.add(piv_id, pivot_trait);
        pmatrix.set(pmatrix_index, child_prob[i]);
    }
}

void Rfunction::evaluate_partner_peel(PeelMatrixKey& pmatrix_index) {
    
    double tmp = 0.0;
    unsigned partner_id;
    enum phased_trait partner_trait;    
    
    partner_id = peel.get_peelnode(0);
    
    for(unsigned i = 0; i < num_alleles; ++i) {
        partner_trait = static_cast<enum phased_trait>(i);
        
        pmatrix_index.add(partner_id, partner_trait);
        
        tmp += (\
                get_disease_probability(partner_id, partner_trait) * \
                previous_rfunction1->get(pmatrix_index) * \
                (previous_rfunction2 == NULL ? 1.0 : previous_rfunction2->get(pmatrix_index)) \
            );
    }
    
    pmatrix.set(pmatrix_index, tmp);
}

void Rfunction::evaluate_element(
                    PeelMatrixKey& pmatrix_index, 
                    DescentGraph* dg, 
                    unsigned locus_index) {
    
    // XXX could remove this with some inheritance?
    switch(peel.get_type()) {
        
        case CHILD_PEEL :
            evaluate_child_peel(pmatrix_index, dg, locus_index);
            break;
            
        case PARTNER_PEEL :
            evaluate_partner_peel(pmatrix_index);
            break;
        
        case PARENT_PEEL :
            evaluate_parent_peel(pmatrix_index, dg, locus_index);
            break;
        
        default :
            fprintf(stderr, "error: default should never be reached! (%s:%d)\n", __FILE__, __LINE__);
            abort();
    }
}

void Rfunction::evaluate(DescentGraph* dg, unsigned int locus_index) {
    PeelMatrixKey k;
    vector<unsigned> q;
    unsigned ndim = peel.get_cutset_size();
    unsigned tmp;
    
    
    // nothing in the cutset to be enumerated
    if(peel.get_type() == LAST_PEEL) {
        evaluate_last_peel(k);
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
            evaluate_element(k, dg, locus_index);
        }
        
        tmp = q.back() + 1;
        q.pop_back();
        
        if(tmp < num_alleles) {
            q.push_back(tmp);
            tmp = ndim - q.size();
            // fill out rest with zeroes
            for(unsigned i = 0; i < tmp; ++i) {
                q.push_back(0);
            }
        }
    }
}

