// generic functions needed by multiple kernels
// eg: handling common r-function things

#include <float.h>

__device__ int gpu_offsets[] = {
    1 <<  0, 
    1 <<  2, 
    1 <<  4, 
    1 <<  6, 
    1 <<  8, 
    1 << 10, 
    1 << 12,
    1 << 14,
    1 << 16,
    1 << 18,
    1 << 20,
    1 << 22,
    1 << 24,
    1 << 26,
    1 << 28,
    1 << 30
};

__device__ float LOG_ZERO = FLT_MIN;

// I am going to assume that the length of 'assignment' is to the number of 
// pedigree members and that everything that is not assigned is -1
//
__device__ int rfunction_index(struct rfunction* rf, int* assignment, int length) {
    int i, tmp = 0;
    
    for(i = 0; i < (rf->cutset_length - 1); ++i) {
        tmp += (assignment[rf->cutset[i]] * gpu_offsets[i]);
    }
    
    return tmp;
}

// this works out the index for the 'presum_matrix', though with this
// scheme I can always just do:
//      rf->matrix[index - ((0..3) * offset[rf->cutset_length - 1])]
// when I sum over one of the dimensions so long as I ensure that
// the peel_node is always in cutset[cutset_length - 1]
__device__ int rfunction_presum_index(struct rfunction* rf, int* assignment, int length) {
    int i, tmp = 0;
    
    for(i = 0; i < rf->cutset_length; ++i) {
        tmp += (assignment[rf->cutset[i]] * gpu_offsets[i]);
    }
    
    return tmp;
}

// from an index, construct the assignment of genotypes to cutset nodes
// 
__device__ void rfunction_presum_assignment(struct rfunction* rf, int ind, int* assignment, int length) {
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

__device__ void rfunction_assignment(struct rfunction* rf, int ind, int* assignment, int length) {
    int index = ind;
    int i;
    
    for(i = 0; i < length; ++i) {
        assignment[i] = -1;
    }
    
    for(i = (rf->cutset_length - 2); i > -1; --i) {
        assignment[rf->cutset[i]] = index / gpu_offsets[i];
        index %= gpu_offsets[i];
    }
}

__device__ float rfunction_get(struct rfunction* rf, int* assignment, int length) {
    return rf->matrix[rfunction_index(rf, assignment, length)];
}

__device__ int get_trait(int value, int parent) {
    switch(parent) {
        case GPU_MATERNAL_ALLELE:
            return (((value == GPU_TRAIT_AA) || (value == GPU_TRAIT_AB)) ? GPU_TRAIT_A : GPU_TRAIT_B);    
        case GPU_PATERNAL_ALLELE:
            return (((value == GPU_TRAIT_AA) || (value == GPU_TRAIT_BA)) ? GPU_TRAIT_A : GPU_TRAIT_B);            
        default:
            break;
    }
    return -1;
}

__device__ float rfunction_trait_prob(struct gpu_state* state, int id, int value, int locus) {
    struct person* p = GET_PERSON(state, id);
    
    if(PERSON_ISTYPED(p)) {
        switch(PERSON_GENOTYPE(p, locus)) {
            case GPU_GENOTYPE_AB :
                return ((value == GPU_TRAIT_AB) || (value == GPU_TRAIT_BA)) ? 0.5 : 0.0;
            case GPU_GENOTYPE_AA :
                return (value == GPU_TRAIT_AA) ? 1.0 : 0.0;
            case GPU_GENOTYPE_BB :
                return (value == GPU_TRAIT_BB) ? 1.0 : 0.0;
            default :
                break;
        }
    }
    
    if(! PERSON_ISFOUNDER(p))
        return 0.25;
    
    return MAP_PROB(GET_MAP(state), locus, value);
}

__device__ float log_sum(float a, float b) {
    if(a == LOG_ZERO)
        return b;
    
    if(b == LOG_ZERO)
        return a;
    
    return log(exp(b - a) + 1) + a;
}

__device__ float log_product(float a, float b) {
    return ((a == LOG_ZERO) or (b == LOG_ZERO)) ? LOG_ZERO : a + b;
}

