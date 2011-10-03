// CUDA locus sampler
// one thread per cell in the pre-sum matrix
// number of threads actually utilised really depends on which 
// r-function is currently being calculated


__device__ double rfunction_trans_prob(struct gpu_state* state, int locus, int peelnode, 
                           int parent_trait, int child_trait, int parent) {
    
    int trait = get_trait(child_trait, parent);
    int meiosis = 0;
    double tmp = 1.0;
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    struct geneticmap* map = GET_MAP(state);
    
    switch(parent_trait) {
        case GPU_TRAIT_AA:
            return (trait == GPU_TRAIT_A) ? 0.25 : 0.0;
        case GPU_TRAIT_BB:
            return (trait == GPU_TRAIT_B) ? 0.25 : 0.0;
        case GPU_TRAIT_AB:
            meiosis = (trait == GPU_TRAIT_A) ? 0 : 1;
            break;
        case GPU_TRAIT_BA:
            meiosis = (trait == GPU_TRAIT_A) ? 1 : 0;
            break;
    }
    
    if(locus != 0) {
        tmp *= (DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus - 1, parent)) == meiosis ? \
                MAP_INVERSETHETA(map, locus - 1) : \
                MAP_THETA(map, locus - 1));
    }
    
    if(locus != (MAP_LENGTH(map) - 1)) {
        tmp *= (DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus + 1, parent)) == meiosis ? \
                MAP_INVERSETHETA(map, locus) : \
                MAP_THETA(map, locus));
    }
    
    return tmp;
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
    double total;
    int i;
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    struct geneticmap* map = GET_MAP(state);
    
    prob_dist[0] = 1.0;
    prob_dist[1] = 1.0;
    
    if(locus != 0) {
        i = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, locus - 1, parent));
        prob_dist[0] *= (i == 0 ? MAP_INVERSETHETA(map, locus - 1) : MAP_THETA(map, locus - 1));
        prob_dist[1] *= (i == 1 ? MAP_INVERSETHETA(map, locus - 1) : MAP_THETA(map, locus - 1));
    }
    
    // map length is actually the number of recombination fractions (#markers - 1)
    if(locus != (MAP_LENGTH(map) - 1)) {
        i = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, locus + 1, parent));
        prob_dist[0] *= (i == 0 ? MAP_INVERSETHETA(map, locus) : MAP_THETA(map, locus));
        prob_dist[1] *= (i == 1 ? MAP_INVERSETHETA(map, locus) : MAP_THETA(map, locus));
    }
    
    total = prob_dist[0] + prob_dist[1];
    
    prob_dist[0] /= total;
    prob_dist[1] /= total;
    
    return (get_random(state) < prob_dist[0]) ? 0 : 1;
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
            //abort();
    }
    
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
    double r = get_random(state);
    int peelnode = RFUNCTION_PEELNODE(rf);
    int i;
    
    // extract probabilities
    for(i = 0; i < NUM_ALLELES; ++i) {
        assignment[peelnode] = i;
        prob_dist[i] = RFUNCTION_PRESUM_GET(rf, rfunction_presum_index(rf, assignment, state->pedigree_length));
        total += prob_dist[i];        
    }
    
    // normalise
    for(i = 0; i < NUM_ALLELES; ++i) {
        prob_dist[i] /= total;
    }
    
    /*
    if((blockIdx.x == 0) && (threadIdx.x == 0)) {
        printf("node = %d\n", RFUNCTION_PEELNODE(rf));
        for(i = 0; i < NUM_ALLELES; ++i) {
            printf(" prob[%d] = %.3f\n", i, prob_dist[i]);
        }
        printf("\n");
    }
    */
    
    // sample
    total = 0.0;
    for(i = 0; i < NUM_ALLELES; ++i) {
        total += prob_dist[i];
        if(r < total) {
            assignment[peelnode] = i;
            return;
        }
    }
    
    //abort();
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
        rfunction_trans_prob(state, locus, peelnode, mother_value, peelnode_value, GPU_MATERNAL_ALLELE) * \
        rfunction_trans_prob(state, locus, peelnode, father_value, peelnode_value, GPU_PATERNAL_ALLELE) * \
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
    
    for(i = 0; i < rf->cutset_length; ++i) {
        p = GET_PERSON(state, rf->cutset[i]);
        
        if(! PERSON_ISPARENT(p, peelnode))
            continue;
            
        tmp *= (rfunction_trans_prob(state, locus, PERSON_ID(p), assignment[PERSON_MOTHER(p)], assignment[PERSON_ID(p)], GPU_MATERNAL_ALLELE) * \
                rfunction_trans_prob(state, locus, PERSON_ID(p), assignment[PERSON_FATHER(p)], assignment[PERSON_ID(p)], GPU_PATERNAL_ALLELE));
    }
    
    RFUNCTION_PRESUM_SET(rf, ind, tmp);
}

__device__ void rfunction_sum(struct rfunction* rf, int ind) {
    int i;
    int tmp = gpu_offsets[rf->cutset_length - 1];
    
    //printf("s b%d t%d\n", blockIdx.x, threadIdx.x);
    
    RFUNCTION_SET(rf, ind, 0.0);
    
    // unroll?
    for(i = 0; i < 4; ++i) {
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
            //abort();
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
    // assumes only one thread active
    for(i = 0; i < rf->presum_length; ++i) {
        rfunction_evaluate_element(rf, state, locus, i);
    }
    
    for(i = 0; i < rf->matrix_length; ++i) {
        rfunction_sum(rf, i);
    }
    */
}

__device__ void sampler_run(struct gpu_state* state, int locus) {
    int i;
    int assignment[128];
    
    // forward peel
    for(i = 0; i < state->functions_per_locus; ++i) {
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

// number of blocks is half the number of loci
__global__ void lsampler_kernel(struct gpu_state* state, int offset) {
    int locus = (blockIdx.x * 2) + offset;

    sampler_run(state, locus);
/*    
    if(locus != (state->map->map_length - 1)) {
        sampler_run(state, locus + 1);
    }
*/
}

void run_gpu_lsampler_kernel(int numblocks, int numthreads, struct gpu_state* state, int offset) {
    lsampler_kernel<<<numblocks, numthreads>>>(state, offset);
}

