using namespace std;

#include <cmath>

#include "rfunction.h"
#include "peeling.h"
#include "peel_matrix.h"
#include "genotype.h"


void Rfunction::generate_key(PeelMatrixKey& pmatrix_index, vector<unsigned int>& assignments) {
    pmatrix_index.reassign(peel.get_cutset(), assignments);
}

void Rfunction::evaluate_element(PeelMatrixKey& pmatrix_index, PeelMatrix* prev_matrix) {
    double tmp;
    double recombination_prob;
    double disease_prob;
    double old_prob;
    enum phased_genotype pg;
    PeelMatrixKey prev_index(pmatrix_index);

    // given that 'prev_matrix' exists, we need to be able to query it
    // how this is performed depends on the 'type' of peel we are talking
    // about and I am not sure whether this procedure is (or can be) particularly
    // general
    //
    // XXX perhaps the PeelMatrixKey class should have the responsibility of 
    // figuring this out?
    switch(peel.get_type()) {
        case CHILD_PEEL :
            // 1. add pivot later
            // 2. remove parents now
            prev_index.remove(pivot->get_maternalid()); // XXX this work should be cached
            prev_index.remove(pivot->get_paternalid()); // XXX this work should be cached
            break;
            
        case PARTNER_PEEL :
            // 1. add pivot later
            break;
        
        case PARENT_PEEL :  // XXX don't bother with yet    
        case LAST_PEEL :    // XXX never seen here? just a sum, handle in 'Rfunction::evaluate'
        default :
            abort();
    }
    
    
    for(unsigned int i = 0; i < num_alleles; ++i) {
        
        //fprintf(stderr, "evaluate_element, allele=%d\n", i);
        
        pg = static_cast<enum phased_genotype>(i);
        
        // finish making the key for the previous rfunction
        prev_index.add(peel.get_pivot(), pg);

        // (a) look up disease prob for this allele with this person
        // (b) look up the recombination prob for this allele, given the DescentGraph + GeneticMap
        // (c) look up relevant old value from previous R function, if 'prev' is not null
        disease_prob = pivot->get_disease_prob(pg);
        recombination_prob = 1.0; // XXX place holder, see b.
        old_prob = prev_matrix != NULL ? prev_matrix->get(prev_index) : 1.0;
        
        // (d) given assignments in 'pmatrix_index' 
        // (sum probabilities for pivot) x (old value from (i)) x (recombination probs)
        tmp += (disease_prob * recombination_prob * old_prob);
    }

    // (e) assign new likelihood in 'pmatrix'
    pmatrix.set(pmatrix_index, tmp);

    pmatrix_index.print();
    printf(" := %f\n", tmp);
}

// XXX can i tell if these matrix can be used together
//
bool Rfunction::evaluate(PeelMatrix* previous_matrix) {
    PeelMatrixKey k;
    vector<unsigned int> q;
    unsigned int ndim = peel.get_cutset_size();
    unsigned int tmp;
    unsigned int i;
        
    // initialise to the first element of matrix
    for(i = 0; i < ndim; ++i) {
        q.push_back(0);
    }

    // enumerate all elements in ndim-dimenstional matrix
    while(not q.empty()) {
        
        if(q.size() == ndim) {
            generate_key(k, q);
            evaluate_element(k, previous_matrix);
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

