#include <cmath>
#include <deque>

#include "rfunction.h"

// given an assignment this returns the identical value that 
// evaluate calculates as 'index'
unsigned int Rfunction::get_index(deque<unsigned int> element) {
    unsigned int index = 0;

    for(unsigned int i = 0; i < element.size(); ++i) {
        index += (pow(num_alleles, i) * element[i]);
    }

    return index;
}

void Rfunction::evaluate_element(deque<unsigned int> element, unsigned int index) {
    double tmp = 0;
    
    for(unsigned int i = 0; i < num_alleles; ++i) {
        // look up disease prob for this allele with this person
        //dp = pivot->get_disease_prob(static_cast<enum phased_genotype>(i));
        // XXX PLACE HOLDER
        // sum of likelihood of all possibilities (each possibility x recombination)
    }

    rfunc[index] = tmp;
}

void Rfunction::evaluate() {
    deque<unsigned int> q;
    unsigned int ndim = peel.get_cutset_size();
    unsigned int tmp;
    unsigned int i;
    unsigned int index = 0;
        
    // initialise to the first element of matrix
    for(i = 0; i < ndim; ++i) {
        q.push_front(0);
    }

    // enumerate all elements in ndim-dimenstional matrix
    while(not q.empty()) {
        
        if(q.size() == ndim) {
            // a complete address for a cell in matrix, where 'q' is 
            // an assignment of alleles for nodes specified by 'cutset'
            // 'index' is where this maps to in 1d array 'rfunc'
            evaluate_element(q, index);
            index++;
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

