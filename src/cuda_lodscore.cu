// calculate LOD scores
// 
// functions derived from TraitRfunction.cc and Rfunction.cc
// but without using an intermediate matrix, just summing everything on
// the fly
// 
// i cannot have a thread per cell in presum matrix, but instead
// one thread per cell in matrix, where each thread calculates the 
// likelihood of each of the four assignments for the peel node
// 


__device__ double lodscore_trait_prob(struct gpu_state* state, int id, int value) {
    struct person* p = GET_PERSON(state, id);
    
    return PERSON_DISEASEPROB(p, value);
}

__device__ double lodscore_trans_prob(struct gpu_state* state, int locus, int peelnode, int m_allele, int p_allele) {
    double tmp = 0.25; // i am only going to remove this term later...
    //double tmp = 1.0;    
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    //struct geneticmap* map = GET_MAP(state);
    //double half_theta = MAP_HALF_THETA(map, locus);
    //double half_inversetheta = MAP_HALF_INVERSETHETA(map, locus);
    
    tmp *= ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus, GPU_MATERNAL_ALLELE)) == m_allele) ? \
            HALF_INVTHETA : HALF_THETA /*half_inversetheta : half_theta*/);
    tmp *= ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus + 1, GPU_MATERNAL_ALLELE)) == m_allele) ? \
            HALF_INVTHETA : HALF_THETA  /*half_inversetheta : half_theta*/);
    
    tmp *= ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus, GPU_PATERNAL_ALLELE)) == p_allele) ? \
            HALF_INVTHETA : HALF_THETA /*half_inversetheta : half_theta*/);
    tmp *= ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus + 1, GPU_PATERNAL_ALLELE)) == p_allele) ? \
            HALF_INVTHETA : HALF_THETA /*half_inversetheta : half_theta*/);
    
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
    //int peelnode_value;
    int i;
    double tmp;
    //printf("e b%d t%d\n", blockIdx.x, threadIdx.x);
    
    rfunction_assignment(rf, ind, assignment, assignment_length);
    
    peelnode = RFUNCTION_PEELNODE(rf);
    //peelnode_value = assignment[peelnode];
    tmp = 0.0;
    
    for(i = 0; i < NUM_ALLELES; ++i) {
        
        assignment[peelnode] = i;
    
        tmp += (\
            lodscore_trait_prob(state, peelnode, i) * \
            (rf->prev1 == NULL ? 1.0 : rfunction_get(rf->prev1, assignment, assignment_length)) * \
            (rf->prev2 == NULL ? 1.0 : rfunction_get(rf->prev2, assignment, assignment_length)) \
        );
    }
    
    RFUNCTION_SET(rf, ind, tmp);
    
    //if(locus == 1)
    //    printf("x %d %e\n", locus, tmp);
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
    double tmp;
    
    //printf("e b%d t%d\n", blockIdx.x, threadIdx.x);
    
    rfunction_assignment(rf, ind, assignment, assignment_length);
    
    peelnode = RFUNCTION_PEELNODE(rf);
    //peelnode_value = assignment[peelnode];
    p = GET_PERSON(state, peelnode);
    mother_value = assignment[PERSON_MOTHER(p)];
    father_value = assignment[PERSON_FATHER(p)];
    tmp = 0.0;
    
    for(i = 0; i < 2; ++i) {
        for(j = 0; j < 2; ++j) {
            
            peelnode_value = get_phased_trait(mother_value, father_value, i, j);
            
            assignment[peelnode] = peelnode_value;
            
            tmp += (\
                lodscore_trait_prob(state, peelnode, peelnode_value) * \
                lodscore_trans_prob(state, locus, peelnode, i, j) * \
                (rf->prev1 == NULL ? 1.0 : rfunction_get(rf->prev1, assignment, assignment_length)) * \
                (rf->prev2 == NULL ? 1.0 : rfunction_get(rf->prev2, assignment, assignment_length)) \
            );
        }
    }
    
    RFUNCTION_SET(rf, ind, tmp);
    
    //if(locus == 1)
    //    printf("c %d %e\n", locus, tmp);
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
    double tmp;
    double child_tmp;
    double child_prob;
    double total;
    
    //printf("e b%d t%d\n", blockIdx.x, threadIdx.x);
    
    rfunction_assignment(rf, ind, assignment, assignment_length);
    
    peelnode = RFUNCTION_PEELNODE(rf);
    //peelnode_value = assignment[peelnode];
    //p = GET_PERSON(state, peelnode);
    total = 0.0;
    
    for(peelnode_value = 0; peelnode_value < NUM_ALLELES; ++peelnode_value) {
    
        assignment[peelnode] = peelnode_value;
    
        tmp = lodscore_trait_prob(state, peelnode, peelnode_value) * \
              (rf->prev1 == NULL ? 1.0 : rfunction_get(rf->prev1, assignment, assignment_length)) * \
              (rf->prev2 == NULL ? 1.0 : rfunction_get(rf->prev2, assignment, assignment_length));
        
        child_prob = 1.0;
        
        for(pid = 0; pid < (rf->cutset_length - 1); ++pid) {
            p = GET_PERSON(state, rf->cutset[pid]);
            
            if(! PERSON_ISPARENT(p, peelnode))
                continue;
            
            mother_value = assignment[PERSON_MOTHER(p)];
            father_value = assignment[PERSON_FATHER(p)];
            
            child_tmp = 0.0;
            
            for(i = 0; i < 2; ++i) {
                for(j = 0; j < 2; ++j) {
                    if(get_phased_trait(mother_value, father_value, i, j) == assignment[PERSON_ID(p)]) {
                        child_tmp += lodscore_trans_prob(state, locus, PERSON_ID(p), i, j);
                    }
                }
            }
            
            child_prob *= child_tmp;
        }
        
        total += (tmp * child_prob);
    }
    
    RFUNCTION_SET(rf, ind, total);
    
    //if(locus == 1)
    //    printf("p %d %e\n", locus, total);
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
    
    for(i = threadIdx.x; i < rf->matrix_length; i += NUM_THREADS) {
        lodscore_evaluate_element(rf, state, locus, i);
    }
    
    /*
    if(threadIdx.x == 0) {
        for(i = 0; i < rf->matrix_length; ++i) {
            lodscore_evaluate_element(rf, state, locus, i);
        }
    }
    */
    
    __syncthreads();
}

// this can be done inline?
__device__ double descentgraph_recombination_prob(struct gpu_state* state, int locus) {
    //struct geneticmap* map = GET_MAP(state);
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    struct person* p;
    //double theta = log(MAP_THETA(map, locus));
    //double antitheta = log(MAP_INVERSETHETA(map, locus));
    double tmp = 0.0;
	int i, j;
	
    for(i = 0; i < state->pedigree_length; ++i) {
        p = GET_PERSON(state, i);
        
        if(PERSON_ISFOUNDER(p))
            continue;
        
        for(j = 0; j < 2; ++j) { // mother and father
            if(DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, i, locus,   j)) != \
               DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, i, locus+1, j))) { 
                tmp += LOG_THETA /*theta*/ ;
            }
            else {
                tmp += LOG_INVTHETA /*antitheta*/ ;
            }
        }
    }
    
    return tmp;
}

__device__ void lodscore_add(struct gpu_state* state, double likelihood) {
    double recomb = descentgraph_recombination_prob(state, blockIdx.x);
    
    //printf("l = %.4f, r = %.4f\n", likelihood, recomb);
    
    state->lodscores[blockIdx.x] = gpu_log_sum(state->lodscores[blockIdx.x], likelihood - recomb - state->dg->transmission_prob);
}

// number of blocks is number of loci - 1
__global__ void lodscore_kernel(struct gpu_state* state) {
    int locus = blockIdx.x;
    int i;
    struct rfunction* rf;
    struct geneticmap* map = GET_MAP(state);
    
    // populate map cache in shared memory
    if(threadIdx.x == 0) {
        map_cache[0] = MAP_HALF_THETA(map, locus);
        map_cache[1] = MAP_HALF_INVERSETHETA(map, locus);
        
        map_cache[2] = log(MAP_THETA(map, locus));
        map_cache[3] = log(MAP_INVERSETHETA(map, locus));
        
        map_length = map->map_length;
    }
    
    __syncthreads();
    
    // forward peel
    for(i = 0; i < state->functions_per_locus; ++i) {
        rf = GET_RFUNCTION(state, i, locus);
        lodscore_evaluate(rf, state, locus);
    }
    
    // get result from last rfunction
    // calculate a lod score
    if(threadIdx.x == 0) {
        lodscore_add(state, log(RFUNCTION_GET(rf, 0)));
    }
}

__global__ void lodscoreinit_kernel(double* lodscores) {
    if(threadIdx.x == 0) {
        lodscores[blockIdx.x] = _LOG_ZERO;
        //printf("[init] %d %.3f\n", blockIdx.x, lodscores[blockIdx.x]);
    }
}

__global__ void lodscorenormalise_kernel(struct gpu_state* state, int count, double trait_likelihood) {
    if(threadIdx.x == 0) {
        //printf("%.4f\n", state->lodscores[blockIdx.x]);
        double tmp = (state->lodscores[blockIdx.x] - log((double)count) - trait_likelihood) / log(10.0);
        state->lodscores[blockIdx.x] = tmp;
    }
}

__global__ void lodscoreprint_kernel(struct gpu_state* state) {
    int i;
    struct geneticmap* map = GET_MAP(state);
    
    for(i = 0; i < (map->map_length - 1); ++i) {
        printf("%d\n\t%.3f\n", i, state->lodscores[i]);
    }
}

void run_gpu_lodscore_kernel(int numblocks, int numthreads, struct gpu_state* state) {
    lodscore_kernel<<<numblocks, numthreads>>>(state);
}

void run_gpu_lodscoreinit_kernel(int numblocks, double* lodscores) {
    lodscoreinit_kernel<<<numblocks, 1>>>(lodscores);
}

void run_gpu_lodscorenormalise_kernel(int numblocks, struct gpu_state* state, int count, double trait_likelihood) {
    lodscorenormalise_kernel<<<numblocks, 1>>>(state, count, trait_likelihood);
}

void run_gpu_lodscoreprint_kernel(struct gpu_state* state) {
    lodscoreprint_kernel<<<1, 1>>>(state);
}

void setup_lodscore_kernel() {
    cudaFuncSetCacheConfig(lodscore_kernel, cudaFuncCachePreferL1);
}

