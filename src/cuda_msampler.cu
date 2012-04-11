#include <stdio.h>

#include "cuda_common.h"


__device__ int founderallele_add(struct gpu_state* state, uint8_t* neighbour_vector, struct adjacent_node2* adjacency_matrix, int mat_fa, int pat_fa, int g) {
    //struct adjacent_node* tmp;
    struct adjacent_node2* tmp;
    int i;
    int neighbours;
    int fa_count;
    
    if(g == GPU_GENOTYPE_UNTYPED) {
        return 1;
    }
    
    fa_count = state->founderallele_count;
    
    // check to see if it exists
    // if it does, then the edge needs to have the same label
    //neighbours = fag->num_neighbours[mat_fa];
    neighbours = neighbour_vector[mat_fa];
    
    for(i = 0; i < neighbours; ++i) {
        //tmp = &(fag->graph[(mat_fa * state->founderallele_count) + i]);
        tmp = &(adjacency_matrix[(mat_fa * fa_count) + i]);
        
        if(tmp->id == pat_fa) {
            return tmp->label == g;
        }
    }
    
    //tmp = &(fag->graph[(mat_fa * state->founderallele_count) + fag->num_neighbours[mat_fa]]);
    tmp = &(adjacency_matrix[(mat_fa * fa_count) + neighbour_vector[mat_fa]]);
    tmp->id = pat_fa;
    tmp->label = g;
    //fag->num_neighbours[mat_fa]++;
    neighbour_vector[mat_fa]++;
    
    if(mat_fa != pat_fa) {
        //tmp = &(fag->graph[(pat_fa * state->founderallele_count) + fag->num_neighbours[pat_fa]]);
        tmp = &(adjacency_matrix[(pat_fa * fa_count) + neighbour_vector[pat_fa]]);
        tmp->id = mat_fa;
        tmp->label = g;
        //fag->num_neighbours[pat_fa]++;
        neighbour_vector[pat_fa]++;
    }
    
    return 1;
}

__device__ int founderallelegraph_populate(struct gpu_state* state, int locus, uint8_t* neighbour_vector, struct adjacent_node2* adjacency_matrix) {
	struct person* p;
	//struct founderallelegraph* fag = GET_FOUNDERALLELEGRAPH(state, locus);
	struct descentgraph* dg = GET_DESCENTGRAPH(state);
	int g;
	int i;
	//int legal = 1;
	int pid;
	int parent_allele;
	int mat, pat;
	__shared__ uint8_t sh_pool[256 * 8];
	uint8_t* founderalleles = &(sh_pool[(locus % 8) * 256]);
	
	// find founder allele assignments, this is only related to the current 
	// descent graph and not whether people are typed or not
	for(i = 0; i < state->pedigree_length; ++i) {
	    pid = state->fa_sequence[i];
	    p = GET_PERSON(state, pid);
	    
	    if(PERSON_ISFOUNDER(p)) {
	        founderalleles[pid * 2] = pid * 2;
	        founderalleles[(pid * 2) + 1] = (pid * 2) + 1;
	    }
	    else {
	        parent_allele = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, pid, locus, GPU_MATERNAL_ALLELE));
	        founderalleles[pid * 2] = founderalleles[ (PERSON_MOTHER(p) * 2) + parent_allele ];
	        
	        parent_allele = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, pid, locus, GPU_PATERNAL_ALLELE));
	        founderalleles[(pid * 2) + 1] = founderalleles[ (PERSON_FATHER(p) * 2) + parent_allele ];
	    }
	}
	
	// construct the actual graph from the assignments and the genotype
	// information
	for(i = 0; i < state->pedigree_length; ++i) {
	    pid = state->fa_sequence[i];
	    p = GET_PERSON(state, pid);
	    
	    if(PERSON_ISTYPED(p)) {
	        mat = founderalleles[pid * 2];
	        pat = founderalleles[(pid * 2) + 1];
	        g = PERSON_GENOTYPE(p, locus);
	        
	        if(! founderallele_add(state, neighbour_vector, adjacency_matrix, mat, pat, g)) {
                //legal = 0;
                return 0;
            }
        }
	}
	
	//return legal;
	return 1;
}

// for loops in the allele graph, hetero is always a contradiction
__device__ int correct_alleles_loop(int g, int allele1) {
    switch(g) {
        case GPU_GENOTYPE_AA:
            return allele1 == 1;
        case GPU_GENOTYPE_BB:
            return allele1 == 2;
        case GPU_GENOTYPE_AB:
            //return 0;
        case GPU_GENOTYPE_UNTYPED:
            //abort();
            break;
    }
    return 0;
}

__device__ int correct_alleles(int g, int allele1) {
    switch(g) {
        case GPU_GENOTYPE_AA:
            return allele1 == 1;
        case GPU_GENOTYPE_BB:
            return allele1 == 2;
        case GPU_GENOTYPE_AB:
            return 1;
        case GPU_GENOTYPE_UNTYPED:
            //abort();
            break;
    }
    return 0;
}

__device__ int correct_alleles(int g, int allele1, int allele2) {
    switch(g) {
        case GPU_GENOTYPE_AA:
            return (allele1 == 1) && (allele2 == 1);
        case GPU_GENOTYPE_BB:
            return (allele1 == 2) && (allele2 == 2);
        case GPU_GENOTYPE_AB:
            return ((allele1 == 1) && (allele2 == 2)) || \
                   ((allele1 == 2) && (allele2 == 1));
        case GPU_GENOTYPE_UNTYPED:
            //abort();
            break;
    }
    return 0;
}


__device__ int legal(struct gpu_state* state, int locus, uint8_t* component, int clength, uint8_t* q, int qlength, uint8_t* neighbour_vector, struct adjacent_node2* adjacency_matrix) {
    int node = component[qlength-1];
    int allele = q[qlength-1];
    int i, j;
    int fa_count = state->founderallele_count;
    
    //struct founderallelegraph* fag = GET_FOUNDERALLELEGRAPH(state, locus);
    //struct adjacent_node* adj;
    struct adjacent_node2* adj;

    //for(i = 0; i < fag->num_neighbours[node]; ++i) {
    for(i = 0; i < neighbour_vector[node]; ++i) {
        //adj = &(fag->graph[(node * state->founderallele_count) + i]);
        adj = &(adjacency_matrix[(node * fa_count) + i]);
        
        // if there is a loop
        if(adj->id == node) {
            if(! correct_alleles_loop(adj->label, allele)) {
                return 0;
            }
            continue;
        }
        
        // test if compatible with label
        if(! correct_alleles(adj->label, allele)) {
            return 0;
        }
        
        // find offset of adjacent node in assignment vector
        for(j = 0; j < clength; ++j) {
            if(component[j] == adj->id) {
                break;
            }
        }
        
        // error if not found
        if(j == clength) {
            // XXX abort on gpu?
            //abort();
            //printf("error in legal()\n");
            __trap();
            //return 0;
        }

        // if not assigned yet, then ignore
        if(j > (qlength-1))
            continue;

        // if assigned, then test legality
        if(! correct_alleles(adj->label, allele, q[j])) {
            return 0;
        }

    }
    
    return 1;
}

__device__ double component_likelihood(struct gpu_state* state, int locus, uint8_t* q, int length) {
    struct geneticmap* map = GET_MAP(state);
    double minor = MAP_MINOR(map, locus);
    double major = MAP_MAJOR(map, locus);
    double tmp = 1.0;
    int i;    
    
    for(i = 0; i < length; ++i) {
        tmp *= ((q[i] == 1) ? major : minor);
    }
    
    return tmp;
}

__device__ double founderallelegraph_enumerate(struct gpu_state* state, int locus, uint8_t* component, int cindex, uint8_t* neighbour_vector, struct adjacent_node2* adjacency_matrix) {
    __shared__ uint8_t sh_pool[128 * 8];
    uint8_t* q = &(sh_pool[(locus % 8) * 128]);
    int qindex = 0;
    int skip = 0;
    double prob = 0.0;
/*
    int i;
    // <debug>
    printf("* ");
    for(i = 0; i < cindex; ++i)
        printf("%d ", component[i]);
    printf("\n");
    // </debug>
*/

    if(cindex == 1) {
        return 1.0;
    }

    while(1) {
        
        while(qindex != cindex) {
            q[qindex++] = 1; // push
            if(! legal(state, locus, component, cindex, q, qindex, neighbour_vector, adjacency_matrix)) {
                skip = 1;
                break;
            }
        }
        
        if(!skip) {
            prob += component_likelihood(state, locus, q, qindex);
        }
        
        while(1) {
            
            while((qindex != 0) and (q[qindex-1] == 2)) { // not empty and last value is 2
                qindex--; // pop
            }
        
            if(qindex == 0) {
                goto no_more_assignments;
            }
        
            q[qindex-1] = 2; // set last value to 2
            
            if(legal(state, locus, component, cindex, q, qindex, neighbour_vector, adjacency_matrix)) {
                skip = 0;
                break;
            }
        }
    }
       
no_more_assignments:
    
    return prob;
}

__device__ double founderallelegraph_likelihood(struct gpu_state* state, int locus, uint8_t* neighbour_vector, struct adjacent_node2* adjacency_matrix) {
    __shared__ uint8_t sh_pool[128 * 3 * 8];
    uint8_t* q         = &(sh_pool[((locus % 8) * 384) +   0]);
    uint8_t* component = &(sh_pool[((locus % 8) * 384) + 128]);
    uint8_t* visited   = &(sh_pool[((locus % 8) * 384) + 256]);
    
    int qindex = 0;
    int cindex;
    int total = state->founderallele_count;
    
    int i;
    int tmp;
    int fa_count;
    
    double tmp_prob;
    double prob = 0.0;
    
    //struct founderallelegraph* fag = GET_FOUNDERALLELEGRAPH(state, locus);
    //struct adjacent_node* tmp2;
    struct adjacent_node2* tmp2;
    
    
    fa_count = state->founderallele_count;
    
    #pragma unroll
    for(i = 0; i < 128; ++i) {
        visited[i] = _WHITE;
    }
    
    do {
        cindex = 0;
        
        // find start point
        for(i = 0; i < fa_count; ++i) {
            if(visited[i] == _WHITE) {
                visited[i] = _GREY;
                q[qindex++] = i;
                break;
            }
        }
        
        while(qindex != 0) {
            tmp = q[--qindex];
            
            //for(i = 0; i < fag->num_neighbours[tmp]; ++i) {
            for(i = 0; i < neighbour_vector[tmp]; ++i) {
                //tmp2 = &(fag->graph[(tmp * fa_count) + i]);
                tmp2 = &(adjacency_matrix[(tmp * fa_count) + i]);
                
                if(visited[tmp2->id] == _WHITE) {
                    visited[tmp2->id] = _GREY;
                    q[qindex++] = tmp2->id;
                }
            }
            
            visited[tmp] = _BLACK;
            total--;
            
            component[cindex++] = tmp;
        }
        
        tmp_prob = founderallelegraph_enumerate(state, locus, component, cindex, neighbour_vector, adjacency_matrix);
		tmp_prob = ((tmp_prob == 0.0) ? _LOG_ZERO : log(tmp_prob));
		
		prob = gpu_log_product(tmp_prob, prob);
    
    } while(total != 0);
    
    return prob;
}

__device__ int founderallele_sample(struct gpu_state* state, struct founderallelegraph* fag) {
    double total;
    
    if((fag->prob[0] == _LOG_ZERO) && (fag->prob[1] == _LOG_ZERO)) {
        //abort();
        //printf("error in founderallele_sample\n");
        //__trap();
    }
    
    if(fag->prob[0] == _LOG_ZERO)
        return 1;
        
    if(fag->prob[1] == _LOG_ZERO)
        return 0;
    
    total = gpu_log_sum(fag->prob[0], fag->prob[1]);
    
    //return (rand() / double(RAND_MAX)) < (fag->prob[0] / total) ? 0 : 1;
    return log(get_random(state)) < (fag->prob[0] - total) ? 0 : 1;
}

__device__ int founderallele_sample2(struct gpu_state* state, float prob0, float prob1) {
    float total;
    
    if((prob0 == _LOG_ZERO) && (prob1 == _LOG_ZERO)) {
        //abort();
        //printf("error in founderallele_sample\n");
        __trap();
    }
    
    if(prob0 == _LOG_ZERO)
        return 1;
        
    if(prob1 == _LOG_ZERO)
        return 0;
    
    total = gpu_log_sum(prob0, prob1);
    
    return logf(get_random(state)) < (prob0 - total) ? 0 : 1;
}

__device__ double founderallele_run(struct gpu_state* state, int locus, int personid, int allele, int value) {
    //struct founderallelegraph* fag = GET_FOUNDERALLELEGRAPH(state, locus);
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    int i;
    int tmp;
    double prob = _LOG_ZERO;
    
    struct adjacent_node2* adjacency_matrix;
    uint8_t* neighbour_vector;
    
    i = state->founderallele_count;
    tmp = (locus % 8) * (((i * i)*2) + i);
    
    neighbour_vector = (uint8_t*) &extern_pool[tmp];
    adjacency_matrix = (struct adjacent_node2*) &extern_pool[tmp + i];
    
    
    // save the previous value
    i = DESCENTGRAPH_OFFSET(dg, personid, locus, allele);
    tmp = DESCENTGRAPH_GET(dg, i);
    DESCENTGRAPH_SET(dg, i, value);
    
    for(i = 0; i < state->founderallele_count; ++i) {
        //fag->num_neighbours[i] = 0;
        neighbour_vector[i] = 0;
    }
    
    if(founderallelegraph_populate(state, locus, neighbour_vector, adjacency_matrix)) {
        prob = founderallelegraph_likelihood(state, locus, neighbour_vector, adjacency_matrix);
    }
    
    //printf("%f\n", prob);
    
    // restore previous value
    i = DESCENTGRAPH_OFFSET(dg, personid, locus, allele);
    DESCENTGRAPH_SET(dg, i, tmp);
    
    return prob;
}

__global__ void msampler_likelihood_kernel(struct gpu_state* state, int meiosis) {
    //int locus = (blockIdx.x * blockDim.x) + threadIdx.x;
    int locus = ((blockIdx.x * 256) + threadIdx.x) / 32;
    int personid = (state->founderallele_count / 2) + (meiosis / 2);
    int allele = meiosis % 2;
    struct founderallelegraph* fag;
    struct geneticmap* map = GET_MAP(state);
    
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    int tmp, tmp2, tmp3;
    
    if(locus < map->map_length) {
        fag = GET_FOUNDERALLELEGRAPH(state, locus);
        
        /*
        fag->prob[0] = founderallele_run(state, locus, personid, allele, 0);
        fag->prob[1] = founderallele_run(state, locus, personid, allele, 1);
        */
        
        if((meiosis == 0) && (allele == 0)) {
            fag->prob[0] = founderallele_run(state, locus, personid, allele, 0);
            fag->prob[1] = founderallele_run(state, locus, personid, allele, 1);
        }
        else {
            tmp = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, locus, allele));
            
            tmp2 = DESCENTGRAPH_GET(dg, 
                    DESCENTGRAPH_OFFSET(dg, (state->founderallele_count / 2) + ((meiosis - 1) / 2), 
                                        locus, (meiosis - 1) % 2));
            
            tmp3 = 1 - tmp; // XXX if i just do this inline i get a invalid write of 8 bytes!
            fag->prob[tmp] = fag->prob[tmp2];
            fag->prob[tmp3] = founderallele_run(state, locus, personid, allele, tmp3);
        }
        
    }
}

__global__ void msampler_sampling_kernel(struct gpu_state* state, int meiosis) {
    int personid = (state->founderallele_count / 2) + (meiosis / 2);
    int allele = meiosis % 2;
    struct geneticmap* map = GET_MAP(state);
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    int i, j;
    
    __shared__ float sh_theta[1024];
    __shared__ float sh_inversetheta[1024];
    __shared__ float sh_matrix[1024][2];
    __shared__ int sh_descentgraph[1024][2];
    
    
    map_length = map->map_length;
    
    // we just have one block for now, 512 threads
    for(i = threadIdx.x; i < map_length; i += 256) {
        sh_descentgraph[i][0] = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, i, GPU_MATERNAL_ALLELE));
        sh_descentgraph[i][1] = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, i, GPU_PATERNAL_ALLELE));
        
        sh_matrix[i][0] = state->graphs[i].prob[0];
        sh_matrix[i][1] = state->graphs[i].prob[1];
        
        if(i < (map_length - 1)) {
            sh_theta[i] = logf((float) MAP_THETA(map, i));
            sh_inversetheta[i] = logf((float) MAP_INVERSETHETA(map, i));
        }
    }
    
    __syncthreads();
    
    
    if(threadIdx.x == 0) {
        // forward
        for(i = 1; i < map_length; ++i) {
            /*
            for(j = 0; j < 2; ++j) {
                sh_matrix[i][j] = gpu_log_product(sh_matrix[i][j], \
                                    gpu_log_sum( \
                                        gpu_log_product(sh_matrix[i-1][j],   sh_theta[i-1]), \
                                        gpu_log_product(sh_matrix[i-1][1-j], sh_inversetheta[i-1]) \
                                    ) \
                                  );
            }
            */
            
            // pragma unroll does not work
            sh_matrix[i][0] = gpu_log_product(sh_matrix[i][0], \
                                    gpu_log_sum( \
                                        gpu_log_product(sh_matrix[i-1][0], sh_theta[i-1]), \
                                        gpu_log_product(sh_matrix[i-1][1], sh_inversetheta[i-1]) \
                                    ) \
                                  );
            sh_matrix[i][1] = gpu_log_product(sh_matrix[i][1], \
                                    gpu_log_sum( \
                                        gpu_log_product(sh_matrix[i-1][1], sh_theta[i-1]), \
                                        gpu_log_product(sh_matrix[i-1][0], sh_inversetheta[i-1]) \
                                    ) \
                                  );
        }
        
        // backward
        i = map_length - 1;
        //DESCENTGRAPH_SET(dg, DESCENTGRAPH_OFFSET(dg, personid, i, allele), founderallele_sample2(state, sh_matrix[i][0], sh_matrix[i][1]));
        sh_descentgraph[i][allele] = founderallele_sample2(state, sh_matrix[i][0], sh_matrix[i][1]);
        
        while(--i >= 0) {
            #pragma unroll
            for(j = 0; j < 2; ++j) {
                //sh_matrix[i][j] = gpu_log_product(sh_matrix[i][j], ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, i+1, allele)) != j) ? sh_theta[i] : sh_inversetheta[i]));
                sh_matrix[i][j] = gpu_log_product(sh_matrix[i][j], ((sh_descentgraph[i+1][allele] != j) ? sh_theta[i] : sh_inversetheta[i]));
                    
            }
            
            //DESCENTGRAPH_SET(dg, DESCENTGRAPH_OFFSET(dg, personid, i, allele), founderallele_sample2(state, sh_matrix[i][0], sh_matrix[i][1]));
            sh_descentgraph[i][allele] = founderallele_sample2(state, sh_matrix[i][0], sh_matrix[i][1]);
        }
    }
    
    __syncthreads();
    
    for(i = threadIdx.x; i < map_length; i += 256) {
        DESCENTGRAPH_SET(dg, DESCENTGRAPH_OFFSET(dg, personid, i, GPU_MATERNAL_ALLELE), sh_descentgraph[i][0]);
        DESCENTGRAPH_SET(dg, DESCENTGRAPH_OFFSET(dg, personid, i, GPU_PATERNAL_ALLELE), sh_descentgraph[i][1]);
    }
}
/*
__global__ void msampler_singlethread_kernel(struct gpu_state* state, int meiosis) {
    int locus;
    int personid = (state->founderallele_count / 2) + (meiosis / 2);
    int allele = meiosis % 2;
    struct founderallelegraph* fag;
    struct geneticmap* map = GET_MAP(state);
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    int i, j;
    
    
    if(threadIdx.x == 0) {        
        // build fags
        for(locus = 0; locus < map->map_length; ++locus) {
            //printf("locus = %d\n", locus);
        
            fag = GET_FOUNDERALLELEGRAPH(state, locus);
            
            for(i = 0; i < 2; ++i)
                fag->prob[i] = founderallele_run(state, locus, personid, allele, i);
            //fag->prob[1] = founderallele_run(state, locus, personid, allele, 1);
        }
        
        // sampling
        
        // forward
        for(i = 1; i < state->map->map_length; ++i) {
            for(j = 0; j < 2; ++j) {
                state->graphs[i].prob[j] = gpu_log_product( \
                                                state->graphs[i].prob[j], \
                                                gpu_log_sum( \
                                                    gpu_log_product(state->graphs[i-1].prob[j],   log(MAP_THETA(map, i-1))), \
                                                    gpu_log_product(state->graphs[i-1].prob[1-j], log(MAP_INVERSETHETA(map, i-1))) \
                                                ) \
                                            );
            }
        }
        
        
        // backward
        i = map->map_length - 1;
        fag = GET_FOUNDERALLELEGRAPH(state, i);
        DESCENTGRAPH_SET(dg, DESCENTGRAPH_OFFSET(dg, personid, i, allele), founderallele_sample(state, fag));
        
        while(--i >= 0) {
            fag = GET_FOUNDERALLELEGRAPH(state, i);
            
            for(j = 0; j < 2; ++j) {
                fag->prob[j] = gpu_log_product(fag->prob[j], ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, i+1, allele)) != j) ? \
                    log(MAP_THETA(map, i)) : \
                    log(MAP_INVERSETHETA(map, i))));
            }
            
            DESCENTGRAPH_SET(dg, DESCENTGRAPH_OFFSET(dg, personid, i, allele), founderallele_sample(state, fag));
        }
    }
}
*/
void run_gpu_msampler_kernel(int numblocks, int numthreads, struct gpu_state* state, int meiosis) {
    //msampler_singlethread_kernel<<<numblocks, numthreads>>>(state, meiosis);
}

void run_gpu_msampler_likelihood_kernel(int numblocks, int numthreads, struct gpu_state* state, int meiosis, size_t shared) {
    msampler_likelihood_kernel<<<numblocks, numthreads, shared>>>(state, meiosis);
}

void run_gpu_msampler_sampling_kernel(struct gpu_state* state, int meiosis) {
    msampler_sampling_kernel<<<1, 256>>>(state, meiosis);
}

void setup_msampler_kernel() {
    cudaFuncSetCacheConfig(msampler_likelihood_kernel, cudaFuncCachePreferShared);
    cudaFuncSetCacheConfig(msampler_sampling_kernel,   cudaFuncCachePreferShared);
}

