// CUDA locus sampler
// one thread per cell in the pre-sum matrix
// number of threads actually utilised really depends on which 
// r-function is currently being calculated


__device__ double rfunction_trans_prob(struct gpu_state* state, int locus, int peelnode, 
                           int parent_trait, int child_trait, int parent) {
    
    int trait = get_trait(child_trait, parent);
    int p = 0;
    double tmp = 1.0;
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    
    // deal with homozygotes first
    if(parent_trait == GPU_TRAIT_AA) {
        tmp = (trait == GPU_TRAIT_A) ? 1.0 : 0.0;
        if(locus != 0)
            tmp *= 0.5;
        if(locus != (map_length - 1))
            tmp *= 0.5;
        return tmp;
    }
    else if(parent_trait == GPU_TRAIT_BB) {
        tmp = (trait == GPU_TRAIT_B) ? 1.0 : 0.0;
        if(locus != 0)
            tmp *= 0.5;
        if(locus != (map_length - 1))
            tmp *= 0.5;
        return tmp;
    }
    
    // heterozygotes are informative, so i can look up
    // the recombination fractions
    if(parent_trait == GPU_TRAIT_AB) {
        p = (trait == GPU_TRAIT_A) ? 0 : 1;
    }
    else {
        p = (trait == GPU_TRAIT_B) ? 0 : 1;
    }
    
    tmp = 0.5;
    
    if(locus != 0) {
        tmp *= ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus-1, parent)) == p) ? \
            INVTHETA_LEFT : THETA_LEFT);
    }
    
    if(locus != (map_length - 1)) {
        tmp *= ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus+1, parent)) == p) ? \
            INVTHETA_RIGHT : THETA_RIGHT);
    }
    
    return tmp;
}

__device__ int trans_cache_index(int mat_trait, int pat_trait, int kid_trait) {
    return (mat_trait * 16) + (pat_trait * 4) + kid_trait;
}

__device__ void transmission_matrix(struct gpu_state* state, int locus, int child_id, double* transmission_matrix) {
    int i;
    
    if(threadIdx.x < 64) {
        int tmp = threadIdx.x;
        int mat_trait = tmp / 16; tmp %= 16;
        int pat_trait = tmp / 4;  tmp %= 4;
        int kid_trait = tmp;
        
        transmission_matrix[threadIdx.x] = \
            rfunction_trans_prob(state, locus, child_id, mat_trait, kid_trait, GPU_MATERNAL_ALLELE) * \
            rfunction_trans_prob(state, locus, child_id, pat_trait, kid_trait, GPU_PATERNAL_ALLELE);
    }
    
    __syncthreads();
    
    if(threadIdx.x < 16) {
        int tmp = threadIdx.x * 4;
        double total = transmission_matrix[tmp]   + transmission_matrix[tmp+1] + 
                       transmission_matrix[tmp+2] + transmission_matrix[tmp+3];
        
        #pragma unroll
        for(i = tmp; i < (tmp + 4); ++i) {
            transmission_matrix[i] /= total;
        }
    }
    
    __syncthreads();
}

__device__ void populate_transmission_matrix(struct rfunction* rf, struct gpu_state* state, int locus) {
    int i;
    
    if(RFUNCTION_TYPE(rf) == GPU_PARENT_PEEL) {
        for(i = 0; i < rf->children_length; ++i) {
            transmission_matrix(state, locus, rf->children[i], &(rf->transmission[i * 64]));
        }
    }
    else if(RFUNCTION_TYPE(rf) == GPU_CHILD_PEEL) {
        transmission_matrix(state, locus, rf->peel_node, rf->transmission);
    }
}

__device__ int sample_hetero_mi(int allele, int trait) {
    if(allele == GPU_TRAIT_A) {
        return (trait == GPU_TRAIT_AB) ? 0 : 1;
    }
    else {
        return (trait == GPU_TRAIT_AB) ? 1 : 0;
    }
}

// find prob of setting mi to 0
// find prob of setting mi to 1
// normalise + sample
__device__ int sample_homo_mi(struct gpu_state* state, int personid, int locus, int parent) {
    double prob_dist[2];
    int i;
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    
    prob_dist[0] = 1.0;
    prob_dist[1] = 1.0;
    
    if(locus != 0) {
        i = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, locus-1, parent));
        prob_dist[0] *= (i == 0 ? INVTHETA_LEFT : THETA_LEFT);
        prob_dist[1] *= (i == 1 ? INVTHETA_LEFT : THETA_LEFT);
    }
    
    if(locus != (map_length - 1)) {
        i = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, locus+1, parent));
        prob_dist[0] *= (i == 0 ? INVTHETA_RIGHT : THETA_RIGHT);
        prob_dist[1] *= (i == 1 ? INVTHETA_RIGHT : THETA_RIGHT);
    }
    
    return (get_random(state) < (prob_dist[0] / (prob_dist[0] + prob_dist[1]))) ? 0 : 1;
}

// if a parent is heterozygous, then there is one choice of meiosis indicator
// if a parent is homozygous, then sample based on meiosis indicators to immediate left and right    
__device__ int sample_mi(struct gpu_state* state, int allele, int trait, int personid, int locus, int parent) {
    switch(trait) {
        case GPU_TRAIT_AB:
        case GPU_TRAIT_BA:
            return sample_hetero_mi(allele, trait);
            
        case GPU_TRAIT_AA:
        case GPU_TRAIT_BB:
            return sample_homo_mi(state, personid, locus, parent);
            
        default:
            break;
    }
    
    //__trap();
    return -1;
}

// sample meiosis indicators
// if a parent is heterozygous, then there is one choice of meiosis indicator
// if a parent is homozygous, then sample based on meiosis indicators to immediate left and right
__device__ void sample_meiosis_indicators(struct gpu_state* state, int* assignment, int locus) {
    struct person* p;
    struct descentgraph* dg;
    int i;
    int mother_trait;
    int father_trait;
    int person_trait;
    int mat_allele;
    int pat_allele;
    int mat_mi;
    int pat_mi;
    
    dg = GET_DESCENTGRAPH(state);
    
    for(i = 0; i < state->pedigree_length; ++i) {
        p = GET_PERSON(state, i);
        
        if(PERSON_ISFOUNDER(p))
            continue;
        
        person_trait = assignment[p->id];
        mother_trait = assignment[p->mother];
        father_trait = assignment[p->father];
        
        mat_allele = (person_trait == GPU_TRAIT_AA) || (person_trait == GPU_TRAIT_AB) ? GPU_TRAIT_A : GPU_TRAIT_B;
        pat_allele = (person_trait == GPU_TRAIT_AA) || (person_trait == GPU_TRAIT_BA) ? GPU_TRAIT_A : GPU_TRAIT_B;
        
        mat_mi = sample_mi(state, mat_allele, mother_trait, i, locus, GPU_MATERNAL_ALLELE);
        pat_mi = sample_mi(state, pat_allele, father_trait, i, locus, GPU_PATERNAL_ALLELE);
        
        DESCENTGRAPH_SET(dg, DESCENTGRAPH_OFFSET(dg, i, locus, GPU_MATERNAL_ALLELE), mat_mi);
        DESCENTGRAPH_SET(dg, DESCENTGRAPH_OFFSET(dg, i, locus, GPU_PATERNAL_ALLELE), pat_mi);
    }
}

__device__ void rfunction_sample(struct rfunction* rf, struct gpu_state* state, int* assignment) {
    double prob_dist[NUM_ALLELES];
    double total = 0.0;
    float r = get_random(state);
    int peelnode = RFUNCTION_PEELNODE(rf);
    int i, last = 0;
    //double tmp;
    
    // extract probabilities
    for(i = 0; i < NUM_ALLELES; ++i) {
        assignment[peelnode] = i;
        prob_dist[i] = RFUNCTION_PRESUM_GET(rf, rfunction_presum_index(rf, assignment, state->pedigree_length));
        total += prob_dist[i];
    }
    
    // normalise
    #pragma unroll
    for(i = 0; i < NUM_ALLELES; ++i) {
        prob_dist[i] /= total;
    }
    
    // sample
    total = 0.0;
    
    for(i = 0; i < NUM_ALLELES; ++i) {
        total += prob_dist[i];
        
        if(r < total) {
            assignment[peelnode] = i;
            return;
        }
        
        if(prob_dist[i] != 0.0) {
            last = i;
        }
    }
    
    assignment[peelnode] = last;
}

__device__ void rfunction_evaluate_partner_peel(struct rfunction* rf, struct gpu_state* state, int locus, int ind) {
    
    int assignment_length = state->pedigree_length;
    int assignment[128];
    int peelnode;
    int peelnode_value;
    
    //printf("e b%d t%d\n", blockIdx.x, threadIdx.x);
    
    rfunction_presum_assignment(rf, ind, assignment, assignment_length);
    
    peelnode = RFUNCTION_PEELNODE(rf);
    peelnode_value = assignment[peelnode];
    
    RFUNCTION_PRESUM_SET(rf, ind, \
        rfunction_trait_prob(state, peelnode, peelnode_value, locus) * \
        (rf->prev1 == NULL ? 1.0 : rfunction_get(rf->prev1, assignment, assignment_length)) * \
        (rf->prev2 == NULL ? 1.0 : rfunction_get(rf->prev2, assignment, assignment_length)));
}

__device__ void rfunction_evaluate_child_peel(struct rfunction* rf, struct gpu_state* state, int locus, int ind) {
    
    struct person* p;
    int assignment_length = state->pedigree_length;
    int assignment[128];
    int peelnode;
    int peelnode_value;
    int mother_value;
    int father_value;
    
    //printf("e b%d t%d\n", blockIdx.x, threadIdx.x);
    
    rfunction_presum_assignment(rf, ind, assignment, assignment_length);
    
    peelnode = RFUNCTION_PEELNODE(rf);
    peelnode_value = assignment[peelnode];
    p = GET_PERSON(state, peelnode);
    mother_value = assignment[PERSON_MOTHER(p)];
    father_value = assignment[PERSON_FATHER(p)];

    RFUNCTION_PRESUM_SET(rf, ind, \
        rfunction_trait_prob(state, peelnode, peelnode_value, locus) * \
        rf->transmission[trans_cache_index(mother_value, father_value, peelnode_value)] * \
        (rf->prev1 == NULL ? 1.0 : rfunction_get(rf->prev1, assignment, assignment_length)) * \
        (rf->prev2 == NULL ? 1.0 : rfunction_get(rf->prev2, assignment, assignment_length)));
}

__device__ void rfunction_evaluate_parent_peel(struct rfunction* rf, struct gpu_state* state, int locus, int ind) {
    
    struct person* p;
    int assignment_length = state->pedigree_length;
    int assignment[128];
    int peelnode;
    int peelnode_value;
    int i;
    double tmp;
    
    //printf("e b%d t%d\n", blockIdx.x, threadIdx.x);
    
    rfunction_presum_assignment(rf, ind, assignment, assignment_length);
    
    peelnode = RFUNCTION_PEELNODE(rf);
    peelnode_value = assignment[peelnode];
    
    tmp = rfunction_trait_prob(state, peelnode, peelnode_value, locus) * \
          (rf->prev1 == NULL ? 1.0 : rfunction_get(rf->prev1, assignment, assignment_length)) * \
          (rf->prev2 == NULL ? 1.0 : rfunction_get(rf->prev2, assignment, assignment_length));
    
    for(i = 0; i < rf->children_length; ++i) {
        p = GET_PERSON(state, rf->children[i]);
        
        tmp *= rf->transmission[(64 * i) + trans_cache_index( \
                                            assignment[PERSON_MOTHER(p)], \
                                            assignment[PERSON_FATHER(p)], \
                                            assignment[PERSON_ID(p)]) \
                                          ];
    }
    
    RFUNCTION_PRESUM_SET(rf, ind, tmp);
}

__device__ void rfunction_sum(struct rfunction* rf, int ind) {
    int i;
    int tmp = gpu_offsets[rf->cutset_length - 1];
    
    //printf("s b%d t%d\n", blockIdx.x, threadIdx.x);
    
    RFUNCTION_SET(rf, ind, 0.0);
    
    #pragma unroll
    for(i = 0; i < NUM_ALLELES; ++i) {
        RFUNCTION_ADD(rf, ind, RFUNCTION_PRESUM_GET(rf, ind + (tmp * i)));
    }
}

__device__ void rfunction_evaluate_element(struct rfunction* rf, struct gpu_state* state, int locus, int ind) {
    
    //printf("e b%d t%d\n", blockIdx.x, threadIdx.x);
    
    switch(RFUNCTION_TYPE(rf)) {
        case GPU_CHILD_PEEL:
            rfunction_evaluate_child_peel(rf, state, locus, ind);
            break;
        case GPU_PARTNER_PEEL:
            rfunction_evaluate_partner_peel(rf, state, locus, ind);
            break;
        case GPU_PARENT_PEEL:
            rfunction_evaluate_parent_peel(rf, state, locus, ind);
            break;
        default:
            //printf("error in rfunction_evaluate_element\n");
            //abort();
            //__trap();
            break;
    }
}

__device__ void rfunction_evaluate(struct rfunction* rf, struct gpu_state* state, int locus) {
    int i;
    
    for(i = threadIdx.x; i < rf->presum_length; i += NUM_THREADS) {
        rfunction_evaluate_element(rf, state, locus, i);
    }
    
    __syncthreads();
    
    for(i = threadIdx.x; i < rf->matrix_length; i += NUM_THREADS) {
        rfunction_sum(rf, i);
    }
    
    __syncthreads();
    
    /*
    if(threadIdx.x < rf->presum_length) {
        rfunction_evaluate_element(rf, state, locus, threadIdx.x);
    }
    
    __syncthreads();
    
    if(threadIdx.x < rf->matrix_length) {
        rfunction_sum(rf, threadIdx.x);
    }
    
    __syncthreads();
    */
}

__global__ void lsampler_kernel(struct gpu_state* state, int window_length, int offset) {
    int i;
    int locus = (blockIdx.x * window_length) + offset;
    int assignment[128];
    
    struct geneticmap* map = GET_MAP(state);
    
    // XXX just in case too many blocks are run
    if(locus > (map->map_length - 1)) {
        return;
    }
    
    // populate map cache in shared memory
    if(threadIdx.x == 0) {
        if(locus != 0) {
            map_cache[0] = MAP_THETA(map, locus - 1);
            map_cache[1] = MAP_INVERSETHETA(map, locus - 1);
        }
        
        if(locus != (map->map_length - 1)) {
            map_cache[2] = MAP_THETA(map, locus);
            map_cache[3] = MAP_INVERSETHETA(map, locus);
        }
        
        map_length = map->map_length;
    }
    
    __syncthreads();
    
    // forward peel
    for(i = 0; i < state->functions_per_locus; ++i) {
        // pre-calculate all transmission probabilities
        populate_transmission_matrix(GET_RFUNCTION(state, i, locus), state, locus);
        
        rfunction_evaluate(GET_RFUNCTION(state, i, locus), state, locus);
    }
    
    // reverse peel, sampling ordered genotypes
    if(threadIdx.x == 0) {
        //printf("sample\n");
        for(i = state->functions_per_locus - 1; i >= 0; --i) {
            rfunction_sample(GET_RFUNCTION(state, i, locus), state, assignment);
        }
        
        sample_meiosis_indicators(state, assignment, locus);
    }
    
    __syncthreads();
}

__global__ void lsampler_onepeel_kernel(struct gpu_state* state, int offset, int function_offset) {
    //int i;
    int locus = (blockIdx.x * 2) + offset;
    //int assignment[128];
    
    struct geneticmap* map = GET_MAP(state);
    
    
    // populate map cache in shared memory
    if(threadIdx.x == 0) {
        if(locus != 0) {
            map_cache[0] = MAP_THETA(map, locus - 1);
            map_cache[1] = MAP_INVERSETHETA(map, locus - 1);
        }
        
        if(locus != (map->map_length - 1)) {
            map_cache[2] = MAP_THETA(map, locus);
            map_cache[3] = MAP_INVERSETHETA(map, locus);
        }
        
        map_length = map->map_length;
    }
    
    __syncthreads();
    
    // forward peel
    //for(i = 0; i < state->functions_per_locus; ++i) {
        populate_transmission_matrix(GET_RFUNCTION(state, function_offset, locus), state, locus);
        rfunction_evaluate(GET_RFUNCTION(state, function_offset, locus), state, locus);
    //}
}    

__global__ void lsampler_sample_kernel(struct gpu_state* state, int offset) {
    int i;
    int locus = (blockIdx.x * 2) + offset;
    int assignment[128];
    struct geneticmap* map = GET_MAP(state);

    // reverse peel, sampling ordered genotypes
    if(threadIdx.x == 0) {
        if(locus != 0) {
            map_cache[0] = MAP_THETA(map, locus - 1);
            map_cache[1] = MAP_INVERSETHETA(map, locus - 1);
        }
        
        if(locus != (map->map_length - 1)) {
            map_cache[2] = MAP_THETA(map, locus);
            map_cache[3] = MAP_INVERSETHETA(map, locus);
        }
        
        map_length = map->map_length;
        
        for(i = state->functions_per_locus - 1; i >= 0; --i) {
            rfunction_sample(GET_RFUNCTION(state, i, locus), state, assignment);
        }
        
        sample_meiosis_indicators(state, assignment, locus);
    }
    
    __syncthreads();
}

void run_gpu_lsampler_kernel(int numblocks, int numthreads, struct gpu_state* state, int window_length, int offset) {
    lsampler_kernel<<<numblocks, numthreads>>>(state, window_length, offset);
}

void run_gpu_lsampler_onepeel_kernel(int numblocks, int numthreads, struct gpu_state* state, int offset, int function_index) {
    lsampler_onepeel_kernel<<<numblocks, numthreads>>>(state, offset, function_index);
}

void run_gpu_lsampler_sample_kernel(int numblocks, int numthreads, struct gpu_state* state, int offset) {
    lsampler_sample_kernel<<<numblocks, numthreads>>>(state, offset);
}

void setup_lsampler_kernel() {
    cudaFuncSetCacheConfig(lsampler_kernel, cudaFuncCachePreferL1);
}

