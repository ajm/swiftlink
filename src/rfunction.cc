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

    pivot = ped->get_by_index(peel.get_pivot());
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

double Rfunction::get_disease_probability(enum phased_trait pt) {
    return pivot->get_disease_prob(pt);
}

// XXX this needs to be somewhere that is not at the same position as 
// a genetic marker
// TODO XXX make this generic so it can do arbitrary points between markers
double Rfunction::get_recombination_probability(
                    DescentGraph* dg, unsigned int locus_index,
                    int maternal_allele, int paternal_allele) {

    double tmp = 1.0;
    double half_recomb_prob;
    
    half_recomb_prob = map->get_theta_halfway(locus_index);
    
    tmp *= dg->get(pivot->get_internalid(), locus_index,   MATERNAL) == maternal_allele ? \
            1.0 - half_recomb_prob : half_recomb_prob;
    tmp *= dg->get(pivot->get_internalid(), locus_index+1, MATERNAL) == maternal_allele ? \
            1.0 - half_recomb_prob : half_recomb_prob;
            
    tmp *= dg->get(pivot->get_internalid(), locus_index,   PATERNAL) == paternal_allele ? \
            1.0 - half_recomb_prob : half_recomb_prob;
    tmp *= dg->get(pivot->get_internalid(), locus_index+1, PATERNAL) == paternal_allele ? \
            1.0 - half_recomb_prob : half_recomb_prob;
    
    return tmp;
}

void Rfunction::evaluate_last_peel(
                    PeelMatrixKey& pmatrix_index, 
                    PeelMatrix* prev_matrix) {  
    
    double tmp = 0.0;
    enum phased_trait pivot_trait;
            
    for(unsigned i = 0; i < num_alleles; ++i) {
        pivot_trait = static_cast<enum phased_trait>(i);
        pmatrix_index.add(pivot->get_internalid(), pivot_trait);
        
        if(prev_matrix and not prev_matrix->is_legal(pmatrix_index)) {
            fprintf(stderr, "key generated is illegal! (%s %d)\n", __FILE__, __LINE__);
            abort();
        }
/*
        printf("\n\ttrait=%d d=%e o=%e\t= %e\n", 
               pivot_trait,
               get_disease_probability(pivot_trait), prev_matrix->get(pmatrix_index), 
               get_disease_probability(pivot_trait) * prev_matrix->get(pmatrix_index));
*/
        tmp += (get_disease_probability(pivot_trait) * prev_matrix->get(pmatrix_index));
    }
        
    pmatrix_index.add(pivot->get_internalid(), (enum phased_trait)0);
    pmatrix.set(pmatrix_index, tmp);
    //printf("\nfinal prob := %e\n", tmp);
}

void Rfunction::evaluate_child_peel(
                    PeelMatrixKey& pmatrix_index, 
                    PeelMatrix* prev_matrix, 
                    DescentGraph* dg,
                    unsigned int locus_index) {

    double tmp = 0.0;
    double recombination_prob;
    double disease_prob;
    double old_prob;
    enum phased_trait piv_trait;
    enum phased_trait mat_trait;
    enum phased_trait pat_trait;
    PeelMatrixKey prev_index(pmatrix_index);
    bool add_pivot = additional.size() == 1;
    
    for(unsigned i = 0; i < missing.size(); ++i) {
        prev_index.remove(missing[i]);
    }
    
    mat_trait = pmatrix_index.get(pivot->get_maternalid());
    pat_trait = pmatrix_index.get(pivot->get_paternalid());
    
    // iterate over all descent graphs to determine child trait 
    // based on parents' traits
    for(int i = 0; i < 2; ++i) {        // maternal
        for(int j = 0; j < 2; ++j) {    // paternal
            
            piv_trait = get_phased_trait(mat_trait, pat_trait, i, j);
            
            if(add_pivot) {
                prev_index.add(pivot->get_internalid(), piv_trait);
            }
            
            if(prev_matrix and not prev_matrix->is_legal(prev_index)) {
                fprintf(stderr, "key generated is illegal! (%s %d)\n", __FILE__, __LINE__);
                abort();
            }
            
            disease_prob        = get_disease_probability(piv_trait);
            recombination_prob  = !dg ? 0.25 : 0.25 * get_recombination_probability(dg, locus_index, i, j);
            old_prob            = prev_matrix != NULL ? prev_matrix->get(prev_index) : 1.0;
            
/*            
            printf("\n\t%d%d trait=%d d=%e r=%e o=%e\t= %e\n", 
                    i, j, piv_trait,
                    disease_prob, recombination_prob, old_prob, 
                    disease_prob * recombination_prob * old_prob);
*/
            tmp += (disease_prob * recombination_prob * old_prob);
        }
    }
    
    pmatrix.set(pmatrix_index, tmp);
    //printf(" := %e\n", tmp);
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
    double disease_prob;
    double recombination_prob;
    double old_prob;
    
    for(unsigned i = 0; i < num_alleles; ++i)
        child_prob[i] = 0.0;
    
    for(unsigned mi = 0; mi < num_alleles; ++mi) {        // mother's genotype
        for(unsigned pi = 0; pi < num_alleles; ++pi) {    // father's genotype
            enum phased_trait m = static_cast<enum phased_trait>(mi);
            enum phased_trait p = static_cast<enum phased_trait>(pi);
            
            for(int i = 0; i < 2; ++i) {                        // maternal inheritance
                for(int j = 0; j < 2; ++j) {                    // paternal inheritance
                    pivot_trait = get_phased_trait(m, p, i, j);
                    
                    prev_index.add(pivot->get_maternalid(), m); // add mother
                    prev_index.add(pivot->get_paternalid(), p); // add father
                    
                    if(prev_matrix and not prev_matrix->is_legal(prev_index)) {
                        fprintf(stderr, "key generated is illegal! (%s %d)\n", __FILE__, __LINE__);
                        abort();
                    }
                    
                    maternal_disease_prob = get_disease_probability(m);
                    paternal_disease_prob = get_disease_probability(p);
                    disease_prob        = get_disease_probability(pivot_trait);
                    recombination_prob  = !dg ? 0.25 : 0.25 * get_recombination_probability(dg, locus_index, i, j);
                    old_prob            = prev_matrix != NULL ? prev_matrix->get(prev_index) : 1.0;
                    
                    child_prob[pivot_trait] += \
                        (maternal_disease_prob * \
                         paternal_disease_prob * \
                         disease_prob * \
                         recombination_prob * \
                         old_prob);
                }
            }
        }
    }
    
    
    for(unsigned i = 0; i < num_alleles; ++i) {
        enum phased_trait pivot_trait = static_cast<enum phased_trait>(i);
        pmatrix_index.add(pivot->get_internalid(), pivot_trait);
        pmatrix.set(pmatrix_index, child_prob[i]);
    }
}

void Rfunction::evaluate_partner_peel(
                    PeelMatrixKey& pmatrix_index, 
                    PeelMatrix* prev_matrix) {
    
    double tmp = 0.0;
    PeelMatrixKey prev_index(pmatrix_index);
    enum phased_trait pivot_trait;
    unsigned pivot_id = pivot->get_internalid();
    
    if(missing.size() != 0) {
        fprintf(stderr, "there should be no missing keys in a 'partner' peel! (%s %d)\n", 
            __FILE__, __LINE__);
        for(unsigned i = 0; i < missing.size(); ++i) {  // XXX temporary
            fprintf(stderr, " %d\n", missing[i]);
        }
        abort();
    }
    
    for(unsigned i = 0; i < num_alleles; ++i) {
        pivot_trait = static_cast<enum phased_trait>(i);
        prev_index.add(pivot_id, pivot_trait);
        
        if(prev_matrix and not prev_matrix->is_legal(prev_index)) {
            fprintf(stderr, "key generated is illegal! (%s %d)\n", __FILE__, __LINE__);
            abort();
        }
/*
        printf("\n\ttrait=%d d=%e o=%e\t= %e\n", 
               pivot_trait,
               get_disease_probability(pivot_trait), prev_matrix->get(prev_index), 
               get_disease_probability(pivot_trait) * prev_matrix->get(prev_index));
*/
        tmp += (get_disease_probability(pivot_trait) * prev_matrix->get(prev_index));
    }
    
    pmatrix.set(pmatrix_index, tmp);
    //printf(" := %e\n", tmp);
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

    
    // ascertain whether previous matrix is compatible with the current matrix
    // given the current r-function being applied to the pedigree
    
    // missing: indices not in the current r-function
    // additional: extra indices in the previous r-function
    missing.clear();
    additional.clear();

    if(not pmatrix.key_intersection(previous_matrix, missing, additional)) {

        //
        // XXX there was obviously some rationale behind all these
        // error conditions, i need to remember what they are and
        // comment accordingly
        //
        
        // missing
        //  CHILD_PEEL if the last peel was child as well, there would be a perfect intersections (this code would not have been reached)
        //  CHILD_PEEL contents of peelset, ie: the child, would be in missing
        // 
        //  PARTNER_PEEL contents of peelset, ie: the partner, would be in missing
        //
        //  PARENT_PEEL contents of peelset, ie: both parents, would be in missing
        // 
        // change to: if missing != peelset
        if(missing.size() > 2) {
            fprintf(stderr, "too many people to remove! (%s %d)\n", __FILE__, __LINE__);
            for(unsigned i = 0; i < missing.size(); ++i) {
                fprintf(stderr, "missing[%d] = %d\n", i, missing[i]);
            }
            abort();
        }

        // addition is of size 1 when CHILD_PEEL, PARTNER_PEEL, it is the person being peeled (child or partner)
        // in PARENT_PEEL it is size 2
        if(additional.size() > 1) { // XXX
            fprintf(stderr, "wrong number of additional keys: %d (%s %d)\n", 
                int(additional.size()), __FILE__, __LINE__);
            for(unsigned i = 0; i < additional.size(); ++i) {
                fprintf(stderr, "additional[%d] = %d\n", i, additional[i]);
            }
            abort();
        }
        
        if((additional.size() == 1) and (additional[0] != pivot->get_internalid())) {
            fprintf(stderr, "additional key is not the pivot! (additional[0] = %d) (%s %d)\n",
                int(additional[0]), __FILE__, __LINE__);
            abort();
        }
        
        // this is the first peel
        if((additional.size() == 0) and not (pivot->isleaf() or pivot->isfounder())) {
            fprintf(stderr, "only leaf nodes do not have additional keys (%s %d)\n",
                __FILE__, __LINE__);
            abort();
        }
    }
    // else if not a child peel abort
    
    
    
    
/*
    printf("\n\nRFUNCTION EVALUATION\n");
    printf("pivot = %d\n", pivot->get_internalid());
    printf("type = %s\n", peel.get_type() == CHILD_PEEL ? "child" : "partner");
    for(unsigned i = 0; i < missing.size(); ++i) {
        printf("missing[%d] = %d\n", i, missing[i]);
    }
    for(unsigned i = 0; i < additional.size(); ++i) {
        printf("additional[%d] = %d\n", i, additional[i]);
    }
*/
    
    
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

