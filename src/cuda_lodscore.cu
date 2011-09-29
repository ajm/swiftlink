#include "cuda_common.h"


__device__ float lodscore_trait_prob(struct gpu_state* state, int id, int value) {
    struct person* p = GET_PERSON(state, id);
    
    return PERSON_DISEASEPROB(p, value);
}

__device__ float lodscore_trans_prob(struct gpu_state* state, int locus, int peelnode, int m_allele, int p_allele) {
    float tmp = 0.25;    
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    struct geneticmap* map = GET_MAP(state);
    float half_theta = MAP_HALF_INVERSETHETA(map, locus);
    float half_inversetheta = MAP_HALF_THETA(map, locus);
    
    tmp *= ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus, GPU_MATERNAL_ALLELE)) == m_allele) ? \
            half_inversetheta : half_theta);
    tmp *= ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus + 1, GPU_MATERNAL_ALLELE)) == m_allele) ? \
            half_inversetheta : half_theta);
    
    tmp *= ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus, GPU_PATERNAL_ALLELE)) == p_allele) ? \
            half_inversetheta : half_theta);
    tmp *= ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus + 1, GPU_PATERNAL_ALLELE)) == p_allele) ? \
            half_inversetheta : half_theta);
    
    return tmp;
}

__device__ int affected_trait(int pt, int allele) {
    
    switch(allele) {
        case 0 :
            return (pt == GPU_TRAIT_BA) || (pt == GPU_TRAIT_BB);
        case 1 :
            return (pt == GPU_TRAIT_AB) || (pt == GPU_TRAIT_BB);
    }
    
    return 0;
}

__device__ int get_phased_trait(int m, int p, int m_allele, int p_allele) {
    int m_affected = affected_trait(m, m_allele);
    int p_affected = affected_trait(p, p_allele);
    int pt;
    
    if(m_affected) {
        pt = p_affected ? GPU_TRAIT_BB : GPU_TRAIT_BA;
    }
    else {
        pt = p_affected ? GPU_TRAIT_AB : GPU_TRAIT_AA;
    }
    
    return pt;
}

__device__ void lodscore_evaluate_partner_peel(struct rfunction* rf, struct gpu_state* state, int locus, int ind) {
    
    int assignment_length = state->pedigree_length;
    int assignment[128];
    int peelnode;
    int peelnode_value;
    int i;
    float tmp;
    //printf("e b%d t%d\n", blockIdx.x, threadIdx.x);
    
    rfunction_assignment(rf, ind, assignment, assignment_length);
    
    peelnode = RFUNCTION_PEELNODE(rf);
    peelnode_value = assignment[peelnode];
    tmp = 0;
    
    for(i = 0; i < NUM_ALLELES; ++i) {
        tmp += \
            rfunction_trait_prob(state, peelnode, peelnode_value, locus) * \
            (rf->prev1 == NULL ? 1.0 : rfunction_get(rf->prev1, assignment, assignment_length)) * \
            (rf->prev2 == NULL ? 1.0 : rfunction_get(rf->prev2, assignment, assignment_length));
    }
    
    RFUNCTION_SET(rf, ind, tmp);
}

__device__ void lodscore_evaluate_child_peel(struct rfunction* rf, struct gpu_state* state, int locus, int ind) {
    
    struct person* p;
    int assignment_length = state->pedigree_length;
    int assignment[128];
    int peelnode;
    int peelnode_value;
    int mother_value;
    int father_value;
    int i, j;
    float tmp;
    
    //printf("e b%d t%d\n", blockIdx.x, threadIdx.x);
    
    rfunction_assignment(rf, ind, assignment, assignment_length);
    
    peelnode = RFUNCTION_PEELNODE(rf);
    //peelnode_value = assignment[peelnode];
    p = GET_PERSON(state, peelnode);
    mother_value = assignment[PERSON_MOTHER(p)];
    father_value = assignment[PERSON_FATHER(p)];
    tmp = 0;
    
    for(i = 0; i < 2; ++i) {
        for(j = 0; j < 2; ++j) {
            
            peelnode_value = get_phased_trait(mother_value, father_value, i, j);
            
            tmp += (\
                lodscore_trait_prob(state, peelnode, peelnode_value) * \
                lodscore_trans_prob(state, locus, peelnode, i, j) * \
                (rf->prev1 == NULL ? 1.0 : rfunction_get(rf->prev1, assignment, assignment_length)) * \
                (rf->prev2 == NULL ? 1.0 : rfunction_get(rf->prev2, assignment, assignment_length)) \
            );
        }
    }
    
    RFUNCTION_SET(rf, ind, tmp);
}

__device__ void lodscore_evaluate_parent_peel(struct rfunction* rf, struct gpu_state* state, int locus, int ind) {
    
    struct person* p;
    int assignment_length = state->pedigree_length;
    int assignment[128];
    int peelnode;
    int peelnode_value;
    int mother_value;
    int father_value;
    int i;
    int j;
    int pid;
    float tmp;
    float child_tmp;
    float child_prob;
    float total;
    
    //printf("e b%d t%d\n", blockIdx.x, threadIdx.x);
    
    rfunction_assignment(rf, ind, assignment, assignment_length);
    
    peelnode = RFUNCTION_PEELNODE(rf);
    //peelnode_value = assignment[peelnode];
    p = GET_PERSON(state, peelnode);
    mother_value = assignment[PERSON_MOTHER(p)];
    father_value = assignment[PERSON_FATHER(p)];
    total = 0.0;
    
    for(peelnode_value = 0; peelnode_value < NUM_ALLELES; ++peelnode_value) {
    
        tmp = rfunction_trait_prob(state, peelnode, peelnode_value, locus) * \
              (rf->prev1 == NULL ? 1.0 : rfunction_get(rf->prev1, assignment, assignment_length)) * \
              (rf->prev2 == NULL ? 1.0 : rfunction_get(rf->prev2, assignment, assignment_length));
        
        child_prob = 1.0;
        
        for(pid = 0; pid < rf->cutset_length; ++pid) {
            p = GET_PERSON(state, rf->cutset[pid]);
            
            if(! PERSON_ISPARENT(p, peelnode))
                continue;
            
            child_tmp = 0.0;
            
            for(i = 0; i < 2; ++i) {
                for(j = 0; j < 2; ++j) {
                    if(get_phased_trait(mother_value, father_value, i, j) == assignment[pid]) {
                        child_tmp += lodscore_trans_prob(state, locus, PERSON_ID(p), i, j);
                    }
                }
            }
            
            child_prob *= child_tmp;
        }
        
        total += (tmp * child_prob);
    }
    
    RFUNCTION_SET(rf, ind, total);
}

__device__ void lodscore_evaluate_element(struct rfunction* rf, struct gpu_state* state, int locus, int ind) {
    
    switch(RFUNCTION_TYPE(rf)) {
        case GPU_CHILD_PEEL:
            lodscore_evaluate_child_peel(rf, state, locus, ind);
            break;
        case GPU_PARTNER_PEEL:
            lodscore_evaluate_partner_peel(rf, state, locus, ind);
            break;
        case GPU_PARENT_PEEL:
            lodscore_evaluate_parent_peel(rf, state, locus, ind);
            break;
        default:
            break;
    }
}

__device__ void lodscore_evaluate(struct rfunction* rf, struct gpu_state* state, int locus) {
    int i;
    
    for(i = threadIdx.x; i < rf->matrix_length; i += 256) {
        lodscore_evaluate_element(rf, state, locus, i);
    }
    
    __syncthreads();
}

// number of blocks is number of loci
__global__ void lodscore_kernel(struct gpu_state* state) {
    int locus = blockIdx.x;
    int i;

    // forward peel
    for(i = 0; i < state->functions_per_locus; ++i) {
        lodscore_evaluate(GET_RFUNCTION(state, i, locus), state, locus);
    }
    
    // get result from last rfunction
    // calculate a lod score
}

void run_gpu_lodscore_kernel(int numblocks, int numthreads, struct gpu_state* state) {
    lodscore_kernel<<<numblocks, numthreads>>>(state);
}

