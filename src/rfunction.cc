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


Rfunction::Rfunction(PeelOperation po, Pedigree* p, GeneticMap* m, unsigned int alleles)
    : pmatrix(po.get_cutset_size(), alleles), 
      peel(po), 
      num_alleles(alleles), 
      map(m),
      ped(p) {
        
    if(alleles != 4)
        abort(); // XXX assumption for now...

    pmatrix.set_keys(peel.get_cutset());
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

void Rfunction::evaluate_last_peel(
                    PeelMatrixKey& pmatrix_index, 
                    PeelMatrix* prev_matrix) {  
    
    double tmp = 0.0;
    enum phased_trait last_trait;
    unsigned last_id = peel.get_peelnode(0);
    
    
    if(peel.get_peelset_size() != 1) {
        fprintf(stderr, "peelset must be size 1 (%s:%d)\n", __FILE__, __LINE__);
        abort();
    }
            
    for(unsigned i = 0; i < num_alleles; ++i) {
        last_trait = static_cast<enum phased_trait>(i);
        pmatrix_index.add(last_id, last_trait);
        
        if(prev_matrix and not prev_matrix->is_legal(pmatrix_index)) {
            fprintf(stderr, "key generated is illegal! (%s %d)\n", __FILE__, __LINE__);
            abort();
        }
        
        tmp += (\
                get_disease_probability(last_id, last_trait) * \
                prev_matrix->get(pmatrix_index) \
            );
    }
    
    pmatrix_index.add(last_id, (enum phased_trait)0);
    pmatrix.set(pmatrix_index, tmp);
}

void Rfunction::evaluate_child_peel(
                    PeelMatrixKey& pmatrix_index, 
                    PeelMatrix* prev_matrix, 
                    DescentGraph* dg,
                    unsigned int locus_index) {

    PeelMatrixKey prev_index(pmatrix_index);
    double recombination_prob;
    double disease_prob;
    double old_prob;
    double tmp = 0.0;
    enum phased_trait kid_trait;
    enum phased_trait mat_trait;
    enum phased_trait pat_trait;
    bool add_child = additional.size() == 1;
    unsigned kid_id = peel.get_peelnode(0);
    Person* kid = ped->get_by_index(kid_id);
    
    
    if(peel.get_peelset_size() != 1) {
        fprintf(stderr, "peelset must be size 1 (%s:%d)\n", __FILE__, __LINE__);
        abort();
    }
    
    for(unsigned i = 0; i < missing.size(); ++i) {
        prev_index.remove(missing[i]);
    }
    
    mat_trait = pmatrix_index.get(kid->get_maternalid());
    pat_trait = pmatrix_index.get(kid->get_paternalid());
    
    // iterate over all descent graphs to determine child trait 
    // based on parents' traits
    for(int i = 0; i < 2; ++i) {        // maternal
        for(int j = 0; j < 2; ++j) {    // paternal
            
            kid_trait = get_phased_trait(mat_trait, pat_trait, i, j);
            
            if(add_child) {
                prev_index.add(kid_id, kid_trait);
            }
            
            if(prev_matrix and not prev_matrix->is_legal(prev_index)) {
                fprintf(stderr, "key generated is illegal! (%s %d)\n", __FILE__, __LINE__);
                abort();
            }
            
            disease_prob        = get_disease_probability(kid_id, kid_trait);
            recombination_prob  = !dg ? 0.25 : 0.25 * get_recombination_probability(dg, locus_index, kid_id, i, j);
            old_prob            = prev_matrix != NULL ? prev_matrix->get(prev_index) : 1.0;
            
            tmp += (disease_prob * recombination_prob * old_prob);
        }
    }
    
    pmatrix.set(pmatrix_index, tmp);
}

void Rfunction::evaluate_parent_peel(
                    PeelMatrixKey& pmatrix_index, 
                    PeelMatrix* prev_matrix, 
                    DescentGraph* dg,
                    unsigned int locus_index) {
    
    // first create a 3d matrix, of mother x father x child
    // child is the pivot
    //
    // do :
    // for gm in mother :
    //   for gf in father :
    //     for all descent graphs :
    //       gc = get_phased_trait
    //       make key
    //       get from prev_matrix value using gm + gf + rest of enumeration
    //       add to what is current in cell (gm, gf, gc)
    // for gc in child
    //   sum everything that has gc in the key
    //
    //
    // ... or:
    // i can do all the above with just a double array of length 4 initialised to zeros...

    // both or neither (in case of prev_matrix == NULL) parents were in the cutset last
    // so always add them and evaluate using a tenary expression
    
    double child_prob[4];
    PeelMatrixKey prev_index(pmatrix_index);
    enum phased_trait pivot_trait;
    double maternal_disease_prob;
    double paternal_disease_prob;
    //double disease_prob;
    double recombination_prob;
    double old_prob;
    
    
    //fprintf(stderr, "evaluate_parent_peel\n");
    
    
    if(peel.get_peelset_size() != 2) {
        fprintf(stderr, "peelset must be size 2 (%s:%d)\n", __FILE__, __LINE__);
        abort();
    }
    
    if(missing.size() != 1) {
        fprintf(stderr, "missing must be size 1, read %d (%s:%d)\n", int(missing.size()), __FILE__, __LINE__);
        abort();
    }
    
    // additional can be length 0,1 or 2
    
    
    for(unsigned i = 0; i < missing.size(); ++i) {
        prev_index.remove(missing[i]);
    }
/*
    for(unsigned i = 0; i < missing.size(); ++i)
        printf("missing: %d\n", int(missing[i]));
    
    for(unsigned i = 0; i < additional.size(); ++i)
        printf("additional: %d\n", int(additional[i]));
*/
    unsigned piv_id = missing[0];
    Person* p = ped->get_by_index(piv_id);
    unsigned mat_id = p->get_maternalid();
    unsigned pat_id = p->get_paternalid();
    
    bool add_mother = find(additional.begin(), additional.end(), mat_id) != additional.end();
    bool add_father = find(additional.begin(), additional.end(), pat_id) != additional.end();
    
    for(unsigned i = 0; i < num_alleles; ++i)
        child_prob[i] = 0.0;
    
    for(unsigned mi = 0; mi < num_alleles; ++mi) {        // mother's genotype
        for(unsigned pi = 0; pi < num_alleles; ++pi) {    // father's genotype
            enum phased_trait m = static_cast<enum phased_trait>(mi);
            enum phased_trait p = static_cast<enum phased_trait>(pi);
            
            for(int i = 0; i < 2; ++i) {                        // maternal inheritance
                for(int j = 0; j < 2; ++j) {                    // paternal inheritance
                    pivot_trait = get_phased_trait(m, p, i, j);
                    
                    if(prev_matrix != NULL) {
                        if(add_mother)
                            prev_index.add(mat_id, m); // add mother
                        if(add_father)
                            prev_index.add(pat_id, p); // add father
                    }
                    
                    //prev_index.print();
                    //printf("\n");
                    
                    if(prev_matrix and not prev_matrix->is_legal(prev_index)) {
                        fprintf(stderr, "key generated is illegal! (%s:%d)\n", __FILE__, __LINE__);
                        abort();
                    }
                    
                    maternal_disease_prob = get_disease_probability(mat_id, m);
                    paternal_disease_prob = get_disease_probability(pat_id, p);
                    //disease_prob        = get_disease_probability(piv_id, pivot_trait); // TODO XXX i think this is used twice! XXX
                    recombination_prob  = !dg ? 0.25 : 0.25 * get_recombination_probability(dg, locus_index, piv_id, i, j);
                    old_prob            = prev_matrix != NULL ? prev_matrix->get(prev_index) : 1.0;
                    
                    child_prob[pivot_trait] += \
                        (maternal_disease_prob * \
                         paternal_disease_prob * \
                         /*disease_prob * \*/
                         recombination_prob * \
                         old_prob);
                }
            }
        }
    }
    
    
    for(unsigned i = 0; i < num_alleles; ++i) {
        enum phased_trait pivot_trait = static_cast<enum phased_trait>(i);
        pmatrix_index.add(piv_id, pivot_trait);
        pmatrix.set(pmatrix_index, child_prob[i]);
    }
}

void Rfunction::evaluate_partner_peel(
                    PeelMatrixKey& pmatrix_index, 
                    PeelMatrix* prev_matrix) {
    
    double tmp = 0.0;
    PeelMatrixKey prev_index(pmatrix_index);
        
    unsigned partner_id;
    enum phased_trait partner_trait;    
    
    if(additional.size() != 1) {
        fprintf(stderr, "PARTNER_PEEL can only peel a single node (tried to peel %d) (%s:%d)\n", int(additional.size()), __FILE__, __LINE__);
        abort();
    }
    
    partner_id = additional[0];
    
    for(unsigned i = 0; i < num_alleles; ++i) {
        partner_trait = static_cast<enum phased_trait>(i);
        
        prev_index.add(partner_id, partner_trait);
        
        if(prev_matrix and not prev_matrix->is_legal(prev_index)) {
            fprintf(stderr, "key generated is illegal! (%s %d)\n", __FILE__, __LINE__);
            abort();
        }
        
        tmp += (\
                get_disease_probability(partner_id, partner_trait) * \
                prev_matrix->get(prev_index)\
            );
    }
    
    pmatrix.set(pmatrix_index, tmp);
}

void Rfunction::evaluate_element(
                    PeelMatrixKey& pmatrix_index, 
                    PeelMatrix* prev_matrix, 
                    DescentGraph* dg, 
                    unsigned int locus_index) {
    
    // given that 'prev_matrix' exists, we need to be able to query it
    // how this is performed depends on the 'type' of peel we are talking
    // about and I am not sure whether this procedure is (or can be) particularly
    // general
    //
    // XXX perhaps the PeelMatrixKey class should have the responsibility of 
    // figuring this out?
    //
    // XXX this could all be sped up with template probably (?)
    switch(peel.get_type()) {
        case CHILD_PEEL :
            evaluate_child_peel(pmatrix_index, prev_matrix, dg, locus_index);
            break;
            
        case PARTNER_PEEL :
            evaluate_partner_peel(pmatrix_index, prev_matrix);
            break;
        
        case PARENT_PEEL :  // XXX don't bother with yet (drop through to abort)
            //fprintf(stderr, "error: PARENT_PEEL is not implemented (%s:%d)\n", __FILE__, __LINE__);
            evaluate_parent_peel(pmatrix_index, prev_matrix, dg, locus_index);
            break;
            
        case LAST_PEEL :    // XXX never seen here
            fprintf(stderr, "error: LAST_PEEL is not dealt with here (%s:%d)\n", __FILE__, __LINE__);
            abort();
            
        default :
            fprintf(stderr, "error: default should never be reached! (%s:%d)\n", __FILE__, __LINE__);
            abort();
    }
}

// XXX can i tell if these matrix can be used together
//
bool Rfunction::evaluate(PeelMatrix* previous_matrix, DescentGraph* dg, unsigned int locus_index) {
    PeelMatrixKey k;
    vector<unsigned int> q;
    unsigned int ndim = peel.get_cutset_size();
    unsigned int tmp;
    unsigned int i;
    
    
    // nothing in the cutset to be enumerated
    if(peel.get_type() == LAST_PEEL) {
        evaluate_last_peel(k, previous_matrix);
        return true;
    }
    
    // missing: indices not in the previous r-function
    // additional: extra indices in the previous r-function
    // this is ascertained by comparing the cutsets (used to index peel_matrices)
    // XXX should this be done in the constructor? can't be there is a different peelmatrix per
    // locus, maybe performed lazily, then never again?
    missing.clear();
    additional.clear();
    
    // XXX check that missing is the same as peel.peelset
    
    pmatrix.key_intersection(previous_matrix, missing, additional);
    
    
    // generate all assignments to iterate through n-dimensional matrix
    
    // initialise to the first element of matrix
    for(i = 0; i < ndim; ++i) {
        q.push_back(0);
    }

    // enumerate all elements in ndim-dimenstional matrix
    while(not q.empty()) {
        
        if(q.size() == ndim) {
            generate_key(k, q);
            //k.print(); //ajm
            evaluate_element(k, previous_matrix, dg, locus_index);
        }
        
        tmp = q.back() + 1;
        q.pop_back();
        
        if(tmp < num_alleles) {
            q.push_back(tmp);
            tmp = ndim - q.size();
            // fill out rest with zeroes
            for(i = 0; i < tmp; ++i) {
                q.push_back(0);
            }
        }
    }
    
    return true;
}

