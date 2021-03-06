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
    //struct person* p = GET_PERSON(state, id);
    
    //return PERSON_DISEASEPROB(p, value);
    
    return disease_prob[value];
}

__device__ double lodscore_trans_prob(struct gpu_state* state, int locus, int peelnode, int m_allele, int p_allele) {
    double tmp = 0.25; // i am only going to remove this term later...
    //double tmp = 1.0;    
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    //struct geneticmap* map = GET_MAP(state);
    //double half_theta = MAP_HALF_THETA(map, locus);
    //double half_inversetheta = MAP_HALF_INVERSETHETA(map, locus);
    
    
    tmp *= ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus, GPU_MATERNAL_ALLELE)) == m_allele) ? \
            INVTHETA_LEFT : THETA_LEFT);
    tmp *= ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus, GPU_PATERNAL_ALLELE)) == p_allele) ? \
            INVTHETA_LEFT : THETA_LEFT);
    
    tmp *= ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus + 1, GPU_MATERNAL_ALLELE)) == m_allele) ? \
            INVTHETA_RIGHT : THETA_RIGHT);
    tmp *= ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, peelnode, locus + 1, GPU_PATERNAL_ALLELE)) == p_allele) ? \
            INVTHETA_RIGHT : THETA_RIGHT);
    
    // this does not work for parent peels
    /*
    tmp *= (descentgraph_left_mat  == m_allele) ? INVTHETA_LEFT  : THETA_LEFT;
    tmp *= (descentgraph_left_pat  == p_allele) ? INVTHETA_LEFT  : THETA_LEFT;
    tmp *= (descentgraph_right_mat == m_allele) ? INVTHETA_RIGHT : THETA_RIGHT;
    tmp *= (descentgraph_right_pat == p_allele) ? INVTHETA_RIGHT : THETA_RIGHT;
    */
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
    int assignment[MAX_PEDIGREE_MEMBERS];
    int peelnode;
    //int peelnode_value;
    int i, j;
    double total;
    double tmp;
    //printf("e b%d t%d\n", blockIdx.x, threadIdx.x);
    
    rfunction_assignment(rf, ind, assignment, assignment_length);
    
    peelnode = RFUNCTION_PEELNODE(rf);
    //peelnode_value = assignment[peelnode];
    total = 0.0;
    
    for(i = 0; i < NUM_ALLELES; ++i) {
        
        assignment[peelnode] = i;
    
        tmp = lodscore_trait_prob(state, peelnode, i);
        
        for(j = 0; j < rf->prev_length; ++j) {
            tmp *= rfunction_get(rf->prev[j], assignment, assignment_length);
        }
        
        total += tmp;
    }
    
    RFUNCTION_SET(rf, ind, total);
}

__device__ void lodscore_evaluate_child_peel(struct rfunction* rf, struct gpu_state* state, int locus, int ind) {
    
    struct person* p;
    int assignment_length = state->pedigree_length;
    int assignment[MAX_PEDIGREE_MEMBERS];
    int peelnode;
    int peelnode_value;
    int mother_value;
    int father_value;
    int i, j, k;
    double tmp;
    double total;
    
    //printf("e b%d t%d\n", blockIdx.x, threadIdx.x);
    
    rfunction_assignment(rf, ind, assignment, assignment_length);
    
    peelnode = RFUNCTION_PEELNODE(rf);
    //peelnode_value = assignment[peelnode];
    p = GET_PERSON(state, peelnode);
    mother_value = assignment[PERSON_MOTHER(p)];
    father_value = assignment[PERSON_FATHER(p)];
    total = 0.0;
    
    for(i = 0; i < 2; ++i) {
        for(j = 0; j < 2; ++j) {
            
            peelnode_value = get_phased_trait(mother_value, father_value, i, j);
            
            assignment[peelnode] = peelnode_value;
            
            tmp = (lodscore_trait_prob(state, peelnode, peelnode_value) * \
                   lodscore_trans_prob(state, locus, peelnode, i, j));
                   
            for(k = 0; k < rf->prev_length; ++k) {
                tmp *= rfunction_get(rf->prev[k], assignment, assignment_length);
            }
            
            total += tmp;
        }
    }
    
    RFUNCTION_SET(rf, ind, total);
}

__device__ void lodscore_evaluate_parent_peel(struct rfunction* rf, struct gpu_state* state, int locus, int ind) {
    
    struct person* p;
    int assignment_length = state->pedigree_length;
    int assignment[MAX_PEDIGREE_MEMBERS];
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
    
        tmp = lodscore_trait_prob(state, peelnode, peelnode_value);
        
        for(i = 0; i < rf->prev_length; ++i) {
            tmp *= rfunction_get(rf->prev[i], assignment, assignment_length);
        }
        
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
    
    for(i = threadIdx.x; i < rf->matrix_length; i += blockDim.x) {
        lodscore_evaluate_element(rf, state, locus, i);
    }
    
    __syncthreads();
    
    
    /*
    // for onepeel kernel
    if(threadIdx.x < rf->matrix_length) {
        lodscore_evaluate_element(rf, state, locus, threadIdx.x);
    }
    
    __syncthreads();
    */
    
    /*
    if(threadIdx.x == 0) {
        for(i = 0; i < rf->matrix_length; ++i) {
            lodscore_evaluate_element(rf, state, locus, i);
        }
    }
    
    __syncthreads();
    */
}

// this can be done inline?
__device__ double descentgraph_recombination_prob(struct gpu_state* state, int locus) {
    struct geneticmap* map = GET_MAP(state);
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    struct person* p;
    //double theta = log(MAP_THETA(map, locus));
    //double antitheta = log(MAP_INVERSETHETA(map, locus));
    double tmp = MAP_THETA(map, locus);
    double log_theta = log(tmp); //log(MAP_THETA(map, locus));
    double log_invtheta = log(1.0 - tmp); //log(MAP_INVERSETHETA(map, locus));
	int i, j;
	
	tmp = 0.0;
	
    for(i = 0; i < state->pedigree_length; ++i) {
        p = GET_PERSON(state, i);
        
        if(PERSON_ISFOUNDER(p))
            continue;
        
        for(j = 0; j < 2; ++j) { // mother and father
            if(DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, i, locus,   j)) != \
               DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, i, locus+1, j))) { 
                tmp += log_theta;
            }
            else {
                tmp += log_invtheta;
            }
        }
    }
    
    return tmp;
}

__device__ void lodscore_add(struct gpu_state* state, int offset, double likelihood) { //, double recomb) {
    double recomb = descentgraph_recombination_prob(state, blockIdx.x);
    int tmp = (state->lodscores_per_marker * blockIdx.x) + offset;
    
    state->lodscores[tmp] = gpu_log_sum_dbl(state->lodscores[tmp], likelihood - recomb - state->dg->transmission_prob);
    
    //printf("GPU prob %f (%f, %f, %f)\n", likelihood - recomb - state->dg->transmission_prob, likelihood, recomb, state->dg->transmission_prob);
    //printf("GPU fin  %f\n", state->lodscores[blockIdx.x]);
}

__device__ void copy_rfunction(struct rfunction* src, struct rfunction* dst) {
    int* from = (int*) src;
    int* to = (int*) dst;
    
    if(threadIdx.x < (sizeof(struct rfunction) / 4)) {
        to[threadIdx.x] = from[threadIdx.x];
    }
    
    __syncthreads();
}

// warning: there is no fail over code yet,
//          though it should be trivial, just have a pointer to main memory if 
//          you run out of cache...
__device__ struct rfunction* cache_rfunction(struct gpu_state* state, int index, int locus, char* cache_mem) {
    struct rfunction* rf = GET_RFUNCTION(state, index, locus);
    struct rfunction* current = (struct rfunction*) &cache_mem[0];
    int i, j;
    char* last_ptr;
    
    copy_rfunction(rf, current);

    current->cutset = (int*)    &current[1];
    current->matrix = (double*) &current->cutset[current->cutset_length];
    current->prev   = (struct rfunction**) &current->matrix[current->matrix_length];
    
    last_ptr = (char*) &current->prev[current->prev_length];
    
    // copy contents of cutset
    if(threadIdx.x < current->cutset_length) {
        current->cutset[threadIdx.x] = rf->cutset[threadIdx.x];
    }
    
    __syncthreads();
    
    
    for(i = 0; i < current->prev_length; ++i) {
        struct rfunction* tmp = rf->prev[i];
        struct rfunction* next = (struct rfunction*) last_ptr;
        
        copy_rfunction(tmp, next);
        
        current->prev[i] = next;
        
        next->cutset = (int*) &next[1];
        next->matrix = (double*) &next->cutset[next->cutset_length];
        
        last_ptr = (char*) &next->matrix[next->matrix_length];
        
        if(threadIdx.x < current->cutset_length) {
            next->cutset[threadIdx.x] = tmp->cutset[threadIdx.x];
        }
        
        __syncthreads();
        
        for(j = threadIdx.x; j < next->matrix_length; j += blockDim.x) {
            next->matrix[j] = tmp->matrix[j];
        }
        
        __syncthreads();
    }
    
    return current;
}

__device__ void write_rfunction(struct gpu_state* state, int index, int locus, struct rfunction* rfsrc) {
    struct rfunction* rfdst = GET_RFUNCTION(state, index, locus);
    int i;
    
    for(i = threadIdx.x; i < rfsrc->matrix_length; i += blockDim.x) {
        rfdst->matrix[i] = rfsrc->matrix[i];
    }
    
    __syncthreads();
}

// number of blocks is number of loci - 1
__global__ void lodscore_kernel(struct gpu_state* state) {
    int locus = blockIdx.x;
    int i;
    struct rfunction* rf;
    struct geneticmap* map = GET_MAP(state);
    //struct descentgraph* dg = GET_DESCENTGRAPH(state);
    struct person* p;
    //double recomb_prob;
    //double theta;
    //double invtheta;
    
#ifdef CUDA_SHAREDMEM_CACHE
    __shared__ char* cache[5120];
#endif
    
    for(int index = 0; index < state->lodscores_per_marker; ++index) {
        // forward peel
        for(i = 0; i < state->functions_per_locus; ++i) {
#ifdef CUDA_SHAREDMEM_CACHE
            rf = cache_rfunction(state, i, locus, cache);
#else
            rf = GET_RFUNCTION(state, i, locus);
#endif
            p = GET_PERSON(state, rf->peel_node);
            
            // populate map cache in shared memory
            if(threadIdx.x == 0) {
                map_cache[0] = MAP_PARTIAL_THETA_LEFT(map, locus, index+1);
                map_cache[1] = 1.0 - map_cache[0];
                
                map_cache[2] = MAP_PARTIAL_THETA_RIGHT(map, locus, index+1);
                map_cache[3] = 1.0 - map_cache[2];
                
                map_length = map->map_length;
                
                //if(index == 0) {
                //    recomb_prob = MAP_THETA(map, locus);
                //    theta = log(recomb_prob);
                //    invtheta = log(1.0 - recomb_prob);
                //    recomb_prob = 0.0;
                //}
            }
            
            // avoid looking up meiosis indicators later...
            //if(threadIdx.x == 0) {
                //descentgraph_left_mat  = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, rf->peel_node, locus,   GPU_MATERNAL_ALLELE));
                //descentgraph_left_pat  = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, rf->peel_node, locus,   GPU_PATERNAL_ALLELE));
                //descentgraph_right_mat = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, rf->peel_node, locus+1, GPU_MATERNAL_ALLELE));
                //descentgraph_right_pat = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, rf->peel_node, locus+1, GPU_PATERNAL_ALLELE));                
            
                //if((index == 0) && ! PERSON_ISFOUNDER(p)) {
                //    recomb_prob += ((descentgraph_left_mat == descentgraph_right_mat) ? invtheta : theta);
                //    recomb_prob += ((descentgraph_left_pat == descentgraph_right_pat) ? invtheta : theta);
                //}
            //}
            
            if(threadIdx.x < 4) {
                disease_prob[threadIdx.x] = PERSON_DISEASEPROB(p, threadIdx.x);
            }
            
            __syncthreads();
            
            lodscore_evaluate(rf, state, locus);

#ifdef CUDA_SHAREDMEM_CACHE
            write_rfunction(state, i, locus, rf);
#endif
        }
        
        // get result from last rfunction
        // calculate a lod score
        if(threadIdx.x == 0) {
            //printf("GPU peel %d %f\n", locus, log(RFUNCTION_GET(rf, 0)));
            lodscore_add(state, index, log(RFUNCTION_GET(rf, 0)));//, recomb_prob);
        }
        
        __syncthreads();
    }
}

__global__ void lodscore_onepeel_kernel(struct gpu_state* state, int function_offset) {
/*
    int locus = blockIdx.x;
    //int i;
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
    //for(i = 0; i < state->functions_per_locus; ++i) {
        rf = GET_RFUNCTION(state, function_offset, locus);
        lodscore_evaluate(rf, state, locus);
    //}
    
    
    // get result from last rfunction
    // calculate a lod score
    if((threadIdx.x == 0) && (function_offset == (state->functions_per_locus - 1))) {
        //printf("GPU peel %d %f\n", locus, log(RFUNCTION_GET(rf, 0)));
        lodscore_add(state, log(RFUNCTION_GET(rf, 0)));
    }
*/
}

__global__ void lodscoreinit_kernel(double* lodscores, int lodscores_per_marker) {
    if(threadIdx.x < lodscores_per_marker) {
        lodscores[(lodscores_per_marker * blockIdx.x) + threadIdx.x] = DBL_LOG_ZERO;
    }
}

__global__ void lodscorenormalise_kernel(struct gpu_state* state, int count, double trait_likelihood) {
    if(threadIdx.x < state->lodscores_per_marker) {
        double tmp = (state->lodscores[(state->lodscores_per_marker * blockIdx.x) + threadIdx.x] - log((double)count) - trait_likelihood) / log(10.0);
        state->lodscores[(state->lodscores_per_marker * blockIdx.x) + threadIdx.x] = tmp;
    }
}

__global__ void lodscoreprint_kernel(struct gpu_state* state) {
    int i;
    struct geneticmap* map = GET_MAP(state);
    
    for(i = 0; i < (map->map_length - 1); ++i) {
        //printf("%d\n\t%.3f\n", i, state->lodscores[i]);
    }
}

void run_gpu_lodscore_kernel(int numblocks, int numthreads, struct gpu_state* state) {
    lodscore_kernel<<<numblocks, numthreads>>>(state);
}

void run_gpu_lodscore_onepeel_kernel(int numblocks, int numthreads, struct gpu_state* state, int function_offset) {
    lodscore_onepeel_kernel<<<numblocks, numthreads>>>(state, function_offset);
}

void run_gpu_lodscoreinit_kernel(int numblocks, int numthreads, double* lodscores, int lodscores_per_marker) {
    lodscoreinit_kernel<<<numblocks, numthreads>>>(lodscores, lodscores_per_marker);
}

void run_gpu_lodscorenormalise_kernel(int numblocks, int numthreads, struct gpu_state* state, int count, double trait_likelihood) {
    lodscorenormalise_kernel<<<numblocks, numthreads>>>(state, count, trait_likelihood);
}

void run_gpu_lodscoreprint_kernel(struct gpu_state* state) {
    lodscoreprint_kernel<<<1, 1>>>(state);
}

void setup_lodscore_kernel() {
    cudaFuncSetCacheConfig(lodscore_kernel, cudaFuncCachePreferL1);
}

