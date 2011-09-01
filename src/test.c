#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NUM_NODES 4
#define NUM_ALLELES 4

enum { 
    TRAIT_AA,
    TRAIT_BA,
    TRAIT_AB,
    TRAIT_BB
};

enum {
    CHILD_PEEL,
    PARTNER_PEEL,
    PARENT_PEEL,
};

int offsets[] = {
    1 <<  0, 
    1 <<  2, 
    1 <<  4, 
    1 <<  6, 
    1 <<  8, 
    1 << 10, 
    1 << 12,
    1 << 14,
    1 << 16,
    1 << 18
};

struct rfunction {
    int peel_node;
    int peel_type;
    
    // THE 'peel_node' MUST ALWAYS BE LOCATED
    // AT cutset[cutset_length - 2], SO IT CAN 
    // BE IGNORED EASILY
    
    int* cutset;                // eg: [1,2,3] (including the 'peel_node')
    int cutset_length;
    
    // how are these two matrices going to be treated?
    // could the index to one be a factor of the other
    // eg: presum_matrix[(x,y,z)] and matrix[(x,y,z) / z]?
    
    float* presum_matrix;
    int presum_length;          // 4 ** (cutset_length - 1)
    
    float* matrix;
    int matrix_length;          // 4 ** cutset_length
    
    struct rfunction* prev1;    // must be NULL if not used
    struct rfunction* prev2;    // must be NULL if not used    
};

#define RFUNCTION_GET(rf_ptr, index)        ((rf)->matrix[(index)])
#define RFUNCTION_SET(rf_ptr, index, value) ((rf)->matrix[(index)] = (value))
#define RFUNCTION_ADD(rf_ptr, index, value) ((rf)->matrix[(index)] += (value))

#define RFUNCTION_PRESUM_GET(rf_ptr, index)         ((rf)->presum_matrix[(index)])
#define RFUNCTION_PRESUM_SET(rf_ptr, index, value)  ((rf)->presum_matrix[(index)] = (value))
#define RFUNCTION_PRESUM_ADD(rf_ptr, index, value)  ((rf)->presum_matrix[(index)] += (value))

void rfunction_zero(struct rfunction* rf) {
    int i;
    
    for(i = 0; i < rf->matrix_length; ++i) {
        rf->matrix[i] = 0.0;
    }
    
    for(i = 0; i < rf->presum_length; ++i) {
        rf->presum_matrix[i] = 0.0;
    }
}

/*
float rfunction_get(struct rfunction* rf, int index) {
    return rf->matrix[index];
}

void rfunction_set(struct rfunction* rf, int index, float val) {
    rf->matrix[index] = val;
}

void rfunction_add(struct rfunction* rf, int index, float val) {
    rf->matrix[index] += val;
}
*/

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
        assignment[rf->cutset[i]] = index / offsets[i];
        index %= offsets[i];
    }
}



// -----

void index2assignment(int ind, int* values, int length) {
    int index = ind;
    int i;
    
    for(i = (length - 1); i > -1; --i) {
        values[i] = index / offsets[i];
        index %= offsets[i];
    }
}

int assignment2index(int* values, int length) {
    int tmp = 0;
    int i;
    
    for(i = 0; i < length; ++i) {
        tmp += (values[i] * offsets[i]);
    }
    
    return tmp;
}

int main(int argc, char** argv) {
    int values[NUM_NODES];
    int i;
    int index;
    
    if(argc != 2) {
        fprintf(stderr, "Usage: %s <index>\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    index = atoi(argv[1]);
    
    index2assignment(index, values, NUM_NODES);
    
    if(index != assignment2index(values, NUM_NODES)) {
        printf("assignment2index() does not work\n");
    }

    printf("%d = [ ", index);
    for(i = 0; i < NUM_NODES; ++i) {
        printf("%d ", values[i]);
    }
    printf("]\n");
    
    return EXIT_SUCCESS;
}

