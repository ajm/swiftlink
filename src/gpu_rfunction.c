#include "gpu_rfunction.h"

// I am going to assume that the length of 'assignment' is to the number of 
// pedigree members and that everything that is not assigned is -1
//
int rfunction_index(struct rfunction* rf, int* assignment, int length) {
    int i, tmp = 0;
    
    for(i = 0; i < (rf->cutset_length - 1); ++i) {
        tmp += (assignment[rf->cutset[i]] * offset[i]);
    }
    
    return tmp;
}

// this works out the index for the 'presum_matrix', though with this
// scheme I can always just do:
//      rf->matrix[index - ((0..3) * offset[rf->cutset_length - 1])]
// when I sum over one of the dimensions so long as I ensure that
// the peel_node is always in cutset[cutset_length - 1]
int rfunction_presum_index(struct rfunction* rf, int* assignment, int length) {
    int i, tmp = 0;
    
    for(i = 0; i < rf->cutset_length; ++i) {
        tmp += (assignment[rf->cutset[i]] * offset[i]);
    }
    
    return tmp;
}

// from an index, construct the assignment of genotypes to cutset nodes
// 
void rfunction_presum_assignment(struct rfunction* rf, int ind, int* assignment, int length) {
    int index = ind;
    int i;
    
    for(i = 0; i < length; ++i) {
        assignment[i] = -1;
    }
    
    for(i = (rf->cutset_length - 1); i > -1; --i) {
        assignment[rf->cutset[i]] = index / gpu_offsets[i];
        index %= gpu_offsets[i];
    }
}

