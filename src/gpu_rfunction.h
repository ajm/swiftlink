#ifndef GPU_RFUNCTION_H_
#define GPU_RFUNCTION_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define NUM_ALLELES 4

enum { 
    GPU_TRAIT_AA,
    GPU_TRAIT_BA,
    GPU_TRAIT_AB,
    GPU_TRAIT_BB
};

enum {
    GPU_GENOTYPE_AA,
    GPU_GENOTYPE_AB,
    GPU_GENOTYPE_BB
};

enum {
    GPU_CHILD_PEEL,
    GPU_PARTNER_PEEL,
    GPU_PARENT_PEEL
};

/*
int gpu_offsets[] = {
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
*/

struct global_state {
    struct rfunction* functions;
    struct person* pedigree;
    
    float* thetas;
    float* antithetas;
    int thetas_length;
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

struct person {
    int id;
    int father;
    int mother;
    
    float prob[4];

    int* genotypes;
    int genotypes_length;
    
    int isfounder;
    int istyped;
};

#define RFUNCTION_GET(rf_ptr, index)        ((rf_ptr)->matrix[(index)])
#define RFUNCTION_SET(rf_ptr, index, value) ((rf_ptr)->matrix[(index)] = (value))
#define RFUNCTION_ADD(rf_ptr, index, value) ((rf_ptr)->matrix[(index)] += (value))

#define RFUNCTION_PRESUM_GET(rf_ptr, index)         ((rf_ptr)->presum_matrix[(index)])
#define RFUNCTION_PRESUM_SET(rf_ptr, index, value)  ((rf_ptr)->presum_matrix[(index)] = (value))
#define RFUNCTION_PRESUM_ADD(rf_ptr, index, value)  ((rf_ptr)->presum_matrix[(index)] += (value))

#define RFUNCTION_PEELNODE(rf_ptr)  ((rf_ptr)->peel_node)
#define RFUNCTION_TYPE(rf_ptr)      ((rf_ptr)->peel_type)

#define THETA(state_ptr, n)     ((state_ptr)->thetas[(n)])
#define ANTITHETA(state_ptr, n) ((state_ptr)->antithetas[(n)])

#define GET_PERSON(state_ptr, n)        ((state_ptr)->pedigree[(n)])

#define PERSON_ISFOUNDER(person_ptr)        ((person_ptr)->isfounder)
#define PERSON_ISTYPED(person_ptr)          ((person_ptr)->istyped)
#define PERSON_DISEASEPROB(person_ptr, n)   ((person_ptr)->prob[(n)])
#define PERSON_GENOTYPE(person_ptr, n)      ((person_ptr)->genotypes[(n)])


int rfunction_index(struct rfunction* rf, int* assignment, int length);
int rfunction_presum_index(struct rfunction* rf, int* assignment, int length);
void rfunction_presum_assignment(struct rfunction* rf, int ind, int* assignment, int length);


#endif
