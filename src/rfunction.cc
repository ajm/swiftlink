using namespace std;

#include <cmath>
#include <deque>

#include "rfunction.h"
#include "peeling.h"
#include "peel_matrix.h"


// assignments is a misnomer, it does not matter which nodes these assignments
// relate to so long as they are used in a consistent manner, ie: this assumes
// that does this r-function, the ordering of the peel operation cutset is 
// constant
void Rfunction::generate_key(PeelMatrixKey& pmatrix_index, deque<unsigned int>& assignments) {
    vector<unsigned int>& cutset = peel.get_cutset();
    
    for(unsigned int i = 0; i < cutset.size(); ++i) {
        pmatrix_index.add(cutset[i], assignments[i]);
    }
}

void Rfunction::evaluate_element(PeelMatrixKey& pmatrix_index) {
    double tmp = 0;
    PeelMatrixKey prev_index(pmatrix_index);

    switch(peel.type) {
        case PEEL_CHILD :
            // 1. add pivot later
            // 2. remove parents now
            prev_index.remove(pivot->maternal_id);
            prev_index.remove(pivot->paternal_id);
            break;
            
        case PEEL_PARTNER :
            // 1. add pivot later
            break;
        
        case PEEL_PARENT :  // XXX don't bother with yet    
        case LAST_PEEL :    // XXX never seen here?
        default :
            abort();
    }
    
    
    for(unsigned int i = 0; i < num_alleles; ++i) {

        // finish making the key for the previous rfunction
        prev_index.add(peel.pivot, i);

        // look up disease prob for this allele with this person
        //dp = pivot->get_disease_prob(static_cast<enum phased_genotype>(i));
        
        // (i)      look up relevant old value, if 'prev' is not null
        // (ii)     given assignments in 'element' (sum probabilities for pivot) x (old value from (i)) x (recombination probs)
        // (iii)    assign new likelihood in 'rfunc'
    }
}

void Rfunction::evaluate() {
    PeelMatrixKey k;
    deque<unsigned int> q;
    unsigned int ndim = peel.get_cutset_size();
    unsigned int tmp;
    unsigned int i;
        
    // initialise to the first element of matrix
    for(i = 0; i < ndim; ++i) {
        q.push_front(0);
    }

    // enumerate all elements in ndim-dimenstional matrix
    while(not q.empty()) {
        
        if(q.size() == ndim) {
            generate_key(k, q);
            evaluate_element(k);
        }
        
        tmp = q.front() + 1;
        q.pop_front();
        
        if(tmp < num_alleles) {
            q.push_front(tmp);
            tmp = ndim - q.size();
            // fill out rest with zeroes
            for(i = 0; i < tmp; ++i) {
                q.push_front(0);
            }
        }
    }
}

