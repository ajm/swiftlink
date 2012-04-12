#include <stdio.h>

#include "cuda_common.h"


// obs is observed genotype (edge), can be anything apart from UNTYPED
// a1 and a2 are alleles (node), can only be HOMOZ_A or HOMOZ_B
__device__ int legal(int obs, int a1, int a2) {
    switch(obs) {
        case GPU_GENOTYPE_AB:
            return a1 != a2;
        case GPU_GENOTYPE_AA:
            return (a1 == GPU_GENOTYPE_AA) && (a2 == GPU_GENOTYPE_AA);
        case GPU_GENOTYPE_BB:
            return (a1 == GPU_GENOTYPE_BB) && (a2 == GPU_GENOTYPE_BB);
        case GPU_GENOTYPE_UNTYPED:
            break;
    }
    
    return 0;
    //__trap();
}

__device__ int get_other_allele(int obs, int a1) {
    switch(obs) {
        case GPU_GENOTYPE_AB:
            return a1 == GPU_GENOTYPE_AA ? GPU_GENOTYPE_BB : GPU_GENOTYPE_AA;
        case GPU_GENOTYPE_AA:
            return a1 == GPU_GENOTYPE_AA ? GPU_GENOTYPE_AA : GPU_GENOTYPE_UNTYPED;
        case GPU_GENOTYPE_BB:
            return a1 == GPU_GENOTYPE_BB ? GPU_GENOTYPE_BB : GPU_GENOTYPE_UNTYPED;
        case GPU_GENOTYPE_UNTYPED:
            break;
    }
    
    return 0;
    //__trap();
}

#define ind(x,y) ((offset * (x)) + (y))
#define get_freq(g) ((g) == GPU_GENOTYPE_AA ? major_freq : minor_freq)

__device__ float founderallelegraph_likelihood(struct gpu_state* state, int locus) {
	struct person* p;
	struct descentgraph* dg = GET_DESCENTGRAPH(state);
	struct geneticmap* map = GET_MAP(state);
	int g;
	int i, j;
	int pid;
	int parent_allele;
	int mat_fa, pat_fa;
    float ret_prob;
    
    
	float minor_freq = MAP_MINOR(map, locus);
	float major_freq = 1.0 - minor_freq;
	
    int tmp0, tmp1;
    tmp0 = (14 * state->founderallele_count) + (2 * state->pedigree_length);
    tmp0 += (tmp0 % 4);
    tmp0 *= (locus % 8);
    
    
    float*   prob               = ( float* ) &extern_pool[tmp0];
    uint8_t* allele_assignment  = (uint8_t*) &prob[2 * state->founderallele_count];
    uint8_t* edge_list          = (uint8_t*) &allele_assignment[2 * state->founderallele_count];
     int8_t* group_membership   = ( int8_t*) &edge_list[2 * state->pedigree_length];
     int8_t* group_fixed        = ( int8_t*) &group_membership[state->founderallele_count];
    uint8_t* group_active       = (uint8_t*) &group_fixed[state->founderallele_count];
    uint8_t* group_size         = (uint8_t*) &group_active[state->founderallele_count];
    
    
    int num_groups = 0;
    int group_index = 0;
    
    int group1, group2;
    int fixed1, fixed2;
    int legal0, legal1, legal2, legal3;
    
    int offset = state->founderallele_count;
    
    
    for(i = 0; i < state->founderallele_count; ++i) {
        group_membership[i] = GPU_DEFAULT_COMPONENT;
        group_fixed[i] = -1;
        group_active[i] = 0;
    }
    
	// find founder allele assignments, this is only related to the current 
	// descent graph and not whether people are typed or not
	for(i = 0; i < state->pedigree_length; ++i) {
	    pid = state->fa_sequence[i];
	    p = GET_PERSON(state, pid);
	    
	    if(PERSON_ISFOUNDER(p)) {
	        edge_list[pid * 2] = pid * 2;
	        edge_list[(pid * 2) + 1] = (pid * 2) + 1;
	    }
	    else {
	        parent_allele = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, pid, locus, GPU_MATERNAL_ALLELE));
	        edge_list[pid * 2] = edge_list[ (PERSON_MOTHER(p) * 2) + parent_allele ];
	        
	        parent_allele = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, pid, locus, GPU_PATERNAL_ALLELE));
	        edge_list[(pid * 2) + 1] = edge_list[ (PERSON_FATHER(p) * 2) + parent_allele ];
	    }
	}
	
    // calculate likelihood, the founder allele graph is only
    // constrained by typed members of the pedigree
    for(i = 0; i < state->pedigree_length; ++i) {
        p = GET_PERSON(state, i);
        if(! PERSON_ISTYPED(p))
            continue;
        
        g = PERSON_GENOTYPE(p, locus);
        
        if(g == GPU_GENOTYPE_UNTYPED)
            continue;
        
        tmp0 = i * 2;
        mat_fa = edge_list[tmp0];
        pat_fa = edge_list[tmp0+1];
        
        // autozygous
        if(mat_fa == pat_fa) {
            if(g == GPU_GENOTYPE_AB) {
                return 0.0;
            }
            
            group1 = group_membership[mat_fa];
            
            if(group1 != GPU_DEFAULT_COMPONENT) {
                // already belongs to a component
                fixed1 = group_fixed[group1];
                if(fixed1 != -1) {
                    if(g != allele_assignment[ind(fixed1, mat_fa)]) {
                        return 0.0;
                    }
                }
                else {
                    if(allele_assignment[ind(0, mat_fa)] == g) {
                        group_fixed[group1] = 0;
                        prob[ind(1, group1)] = 0.0;
                    }
                    else if(allele_assignment[ind(1, mat_fa)] == g) {
                        group_fixed[group1] = 1;
                        prob[ind(0, group1)] = 0.0;
                    }
                    else {
                        return 0.0;
                    }
                }
            }
            else {
                group_membership[mat_fa] = group_index;
                group_fixed[group_index] = 0;
                group_active[group_index] = 1; // true
                group_size[group_index] = 1;
                allele_assignment[ind(0, mat_fa)] = g;
                prob[ind(0, group_index)] = get_freq(g);
                prob[ind(1, group_index)] = 0.0;
                
                ++group_index;
                ++num_groups;
            }
        }
        // not autozygous
        else {
            group1 = group_membership[mat_fa];
            group2 = group_membership[pat_fa];
            
            if(group1 != GPU_DEFAULT_COMPONENT) {
                
                fixed1 = group_fixed[group1];
                
                if(group2 != GPU_DEFAULT_COMPONENT) {
                    // both in the same group
                    if(group1 == group2) {
                        if(fixed1 != -1) {
                            // fixed, check if legit
                            if(! legal(g, allele_assignment[ind(fixed1, mat_fa)], allele_assignment[ind(fixed1, pat_fa)])) {
                                return 0.0;
                            }
                        }
                        else {
                            // not fixed, check if still the case
                            legal0 = legal(g, allele_assignment[ind(0, mat_fa)], allele_assignment[ind(0, pat_fa)]);
                            legal1 = legal(g, allele_assignment[ind(1, mat_fa)], allele_assignment[ind(1, pat_fa)]);
                            
                            if(legal0) {
                                if(! legal1) {
                                    group_fixed[group1] = 0;
                                    prob[ind(1, group1)] = 0.0;
                                }
                            }
                            else {
                                if(legal1) {
                                    group_fixed[group1] = 1;
                                    prob[ind(0, group1)] = 0.0;
                                }
                                else {
                                    return 0.0;
                                }
                            }
                        }
                    }
                    else {
                        // in different groups
                        fixed2 = group_fixed[group2];
                        
                        if(fixed1 != -1) {
                            if(fixed2 != -1) {
                                // fixed, check if legit
                                if(! legal(g, allele_assignment[ind(fixed1, mat_fa)], allele_assignment[ind(fixed2, pat_fa)])) {
                                    return 0.0;
                                }
                            }
                            else {
                                // group1 is fixed, which assignment for group 2 works
                                legal0 = legal(g, allele_assignment[ind(fixed1, mat_fa)], allele_assignment[ind(0, pat_fa)]);
                                legal1 = legal(g, allele_assignment[ind(fixed1, mat_fa)], allele_assignment[ind(1, pat_fa)]);
                                
                                if(legal0)
                                    fixed2 = 0;
                                else if(legal1)
                                    fixed2 = 1;
                                else {
                                    return 0.0;
                                }
                            }
                        }
                        else if(fixed2 != -1) {
                            // group2 is fixed, which assignment for group 1 works
                            legal0 = legal(g, allele_assignment[ind(0, mat_fa)], allele_assignment[ind(fixed2, pat_fa)]);
                            legal1 = legal(g, allele_assignment[ind(1, mat_fa)], allele_assignment[ind(fixed2, pat_fa)]);
                            
                            if(legal0)
                                fixed1 = 0;
                            else if(legal1)
                                fixed1 = 1;
                            else {
                                return 0.0;
                            }
                        }
                        else {
                            // neither group1 nor group2 are fixed
                            legal0 = legal(g, allele_assignment[ind(0, mat_fa)], allele_assignment[ind(0, pat_fa)]);
                            legal1 = legal(g, allele_assignment[ind(1, mat_fa)], allele_assignment[ind(0, pat_fa)]);
                            legal2 = legal(g, allele_assignment[ind(0, mat_fa)], allele_assignment[ind(1, pat_fa)]);
                            legal3 = legal(g, allele_assignment[ind(1, mat_fa)], allele_assignment[ind(1, pat_fa)]);
                            
                            if(! (legal0 || legal1 || legal2 || legal3)) {
                                return 0.0;
                            }
                            
                            if(legal0 && !(legal1 || legal2 || legal3)) {
                                // fixed, allele 0, no swap
                                fixed1 = fixed2 = 0;
                            }
                            else if(legal1 && !(legal0 || legal2 || legal3)) {
                                // fixed, swap assignment in one group
                                fixed1 = 1;
                                fixed2 = 0;
                            }
                            else if(legal2 && !(legal0 || legal1 || legal3)) {
                                // fixed, swap assignment in one group
                                fixed1 = 0;
                                fixed2 = 1;
                            }
                            else if(legal3 && !(legal0 || legal1 || legal2)) {
                                // fixed, allele 1, no swap
                                fixed1 = fixed2 = 1;
                            }
                            else if(legal0 && legal3 && !(legal1 || legal2)){
                                // still unfixed, assignments don't need swapping
                                fixed1 = fixed2 = -1;
                            }
                            else if(legal1 && legal2 && !(legal0 || legal3)){
                                // still unfixed, swap assignment in one group
                                fixed1 = fixed2 = -2;
                            }
                            else {
                                // all true, assignments don't need swapping
                                fixed1 = fixed2 = -1;
                            }
                        }
                        
                        int flip = 1;
                        
                        if(fixed1 != fixed2) {
                            // XXX for now, always keep group1
                            group_fixed[group1] = fixed1;
                            prob[ind(1-fixed1, group1)] = 0.0;
                            prob[ind(1-fixed2, group2)] = 0.0;
                            //combine_components(group1, group2, 1);
                        }
                        else {
                            if(fixed1 == -2) {
                                group_fixed[group1] = -1;
                                //combine_components(group1, group2, 1);
                            }
                            else {
                                group_fixed[group1] = fixed1;
                                if(fixed1 != -1) {
                                    prob[ind(1-fixed1, group1)] = 0.0;
                                    prob[ind(1-fixed2, group2)] = 0.0;
                                }
                                //combine_components(group1, group2, 0);
                                flip = 0;
                            }
                        }
                        
                        
                        // XXX begin combine component code
                        for(j = 0; j < state->founderallele_count; ++j) {
                            if(group_membership[j] == group2) {
                                group_membership[j] = group1;
                                
                                if(flip) {
                                    tmp0 = allele_assignment[ind(0, j)];
                                    allele_assignment[ind(0, j)] = allele_assignment[ind(1, j)];
                                    allele_assignment[ind(1, j)] = tmp0;
                                }
                            }
                        }
                        
                        if(flip) {
                            prob[ind(0, group1)] *= prob[ind(1, group2)];
                            prob[ind(1, group1)] *= prob[ind(0, group2)];
                        }
                        else {
                            prob[ind(0, group1)] *= prob[ind(0, group2)];
                            prob[ind(1, group1)] *= prob[ind(1, group2)];
                        }
                        
                        group_size[group1] += group_size[group2];
                        group_active[group2] = 0; // false
                        // XXX end combine component code
                        
                        
                        --num_groups;
                    }
                }
                else {
                    if(fixed1 != -1) {
                        tmp0 = get_other_allele(g, allele_assignment[ind(fixed1, mat_fa)]);
                        
                        if(tmp0 == GPU_GENOTYPE_UNTYPED) {
                            return 0.0;
                        }
                        else {
                            allele_assignment[ind(fixed1, pat_fa)] = tmp0;
                            prob[ind(fixed1, group1)] *= get_freq(tmp0);
                        }
                    }
                    else {
                        tmp0 = get_other_allele(g, allele_assignment[ind(0, mat_fa)]);
                        tmp1 = get_other_allele(g, allele_assignment[ind(1, mat_fa)]);
                        
                        if(tmp0 != GPU_GENOTYPE_UNTYPED) {
                            if(tmp1 != GPU_GENOTYPE_UNTYPED) {
                                allele_assignment[ind(0, pat_fa)] = tmp0;
                                allele_assignment[ind(1, pat_fa)] = tmp1;
                                prob[ind(0, group1)] *= get_freq(tmp0);
                                prob[ind(1, group1)] *= get_freq(tmp1);
                            }
                            else {
                                allele_assignment[ind(0, pat_fa)] = tmp0;
                                prob[ind(0, group1)] *= get_freq(tmp0);
                                prob[ind(1, group1)] = 0.0;
                                group_fixed[group1] = 0;
                            }
                        }
                        else {
                            if(tmp1 != GPU_GENOTYPE_UNTYPED) {
                                allele_assignment[ind(1, pat_fa)] = tmp1;
                                prob[ind(1, group1)] *= get_freq(tmp1);
                                prob[ind(0, group1)] = 0.0;
                                group_fixed[group1] = 1;
                            }
                            else {
                                return 0.0;
                            }
                        }
                    }
                    
                    group_membership[pat_fa] = group1;
                    group_size[group1] += 1;
                }
            }
            else if(group2 != GPU_DEFAULT_COMPONENT) {
                
                fixed2 = group_fixed[group2];
                
                if(fixed2 != -1) {
                    tmp0 = get_other_allele(g, allele_assignment[ind(fixed2, pat_fa)]);
                    
                    if(tmp0 == GPU_GENOTYPE_UNTYPED) {
                        return 0.0;
                    }
                    else {
                        allele_assignment[ind(fixed2, mat_fa)] = tmp0;
                        prob[ind(fixed2, group2)] *= get_freq(tmp0);
                    }
                }
                else {
                    tmp0 = get_other_allele(g, allele_assignment[ind(0, pat_fa)]);
                    tmp1 = get_other_allele(g, allele_assignment[ind(1, pat_fa)]);
                    
                    if(tmp0 != GPU_GENOTYPE_UNTYPED) {
                        if(tmp1 != GPU_GENOTYPE_UNTYPED) {
                            allele_assignment[ind(0, mat_fa)] = tmp0;
                            allele_assignment[ind(1, mat_fa)] = tmp1;
                            prob[ind(0, group2)] *= get_freq(tmp0);
                            prob[ind(1, group2)] *= get_freq(tmp1);
                        }
                        else {
                            allele_assignment[ind(0, mat_fa)] = tmp0;
                            prob[ind(0, group2)] *= get_freq(tmp0);
                            prob[ind(1, group2)] = 0.0;
                            group_fixed[group2] = 0;
                        }
                    }
                    else {
                        if(tmp1 != GPU_GENOTYPE_UNTYPED) {
                            allele_assignment[ind(1, mat_fa)] = tmp1;
                            prob[ind(1, group2)] *= get_freq(tmp1);
                            prob[ind(0, group2)] = 0.0;
                            group_fixed[group2] = 1;
                        }
                        else {
                            return 0.0;
                        }
                    }
                }
                
                group_membership[mat_fa] = group2;
                group_size[group2] += 1;
            }
            else {
                if(g == GPU_GENOTYPE_AB) {
                    allele_assignment[ind(0, mat_fa)] = allele_assignment[ind(1, pat_fa)] = GPU_GENOTYPE_AA;
                    allele_assignment[ind(1, mat_fa)] = allele_assignment[ind(0, pat_fa)] = GPU_GENOTYPE_BB;
                    prob[ind(0, group_index)] = \
                        prob[ind(1, group_index)] = major_freq * minor_freq;
                    group_fixed[group_index] = -1;
                }
                else {
                    allele_assignment[ind(0, mat_fa)] = allele_assignment[ind(0, pat_fa)] = g;
                    prob[ind(0, group_index)] = pow(get_freq(g), 2);
                    prob[ind(1, group_index)] = 0.0;
                    group_fixed[group_index] = 0;
                }
                
                // neither in a group
                group_membership[mat_fa] = \
                    group_membership[pat_fa] = group_index;
                group_size[group_index] = 2;
                group_active[group_index] = 1; // true
                ++group_index;
                ++num_groups;
            }
        }
    }
    
    
    // calculate likelihood
    ret_prob = 1.0;
    
    for(i = 0; i < group_index; ++i) {
        if(group_active[i] /*&& (group_size[i] > 1)*/) {
            fixed1 = group_fixed[i];
            if(fixed1 != -1) {
                ret_prob *= prob[ind(fixed1, i)];
            }
            else {
                ret_prob *= (prob[ind(0, i)] + prob[ind(1, i)]);
            }
        }
    }
    
    return ret_prob;
}

__device__ int founderallele_sample2(struct gpu_state* state, float prob0, float prob1) {
    
    //printf("sample : %e %e\n", prob0, prob1);
    
    if((prob0 == 0.0) && (prob1 == 0.0)) {
        //printf("error in founderallele_sample\n");
        __trap();
    }
    
    return get_random(state) < (prob0 / (prob0 + prob1)) ? 0 : 1;
}

__device__ double founderallele_run(struct gpu_state* state, int locus, int personid, int allele, int value) {
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    int i, tmp;
    double prob = 0.0;
    
    // save the previous value
    i = DESCENTGRAPH_OFFSET(dg, personid, locus, allele);
    tmp = DESCENTGRAPH_GET(dg, i);
    DESCENTGRAPH_SET(dg, i, value);
    
    prob = founderallelegraph_likelihood(state, locus);
    
    // restore previous value
    i = DESCENTGRAPH_OFFSET(dg, personid, locus, allele);
    DESCENTGRAPH_SET(dg, i, tmp);
    
    return prob;
}

__global__ void msampler_reset_kernel(struct gpu_state* state, int meiosis) {
    int locus = ((blockIdx.x * 256) + threadIdx.x) / 32;
    int personid = (state->founderallele_count / 2) + (meiosis / 2);
    int allele = meiosis % 2;
    struct geneticmap* map = GET_MAP(state);
    
    //struct descentgraph* dg = GET_DESCENTGRAPH(state);
    //int tmp, tmp2, tmp3;
    int index = locus * 2;
    
    if(locus < map->map_length) {        
        state->raw_matrix[index]     = founderallele_run(state, locus, personid, allele, 0);
        state->raw_matrix[index + 1] = founderallele_run(state, locus, personid, allele, 1);
    }
}

__global__ void msampler_likelihood_kernel(struct gpu_state* state, int meiosis, int last_meiosis) {
    int locus = ((blockIdx.x * 256) + threadIdx.x) / 32;
    int personid = (state->founderallele_count / 2) + (meiosis / 2);
    int allele = meiosis % 2;
    
    int last_personid = (state->founderallele_count / 2) + (last_meiosis / 2);
    int last_allele = last_meiosis % 2;
    
    
    struct geneticmap* map = GET_MAP(state);
    
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    int tmp, tmp2, tmp3;
    int index = locus * 2;
    
    if(locus < map->map_length) {
        /*
        if((meiosis == 0) && (allele == 0)) {
            state->raw_matrix[index]     = founderallele_run(state, locus, personid, allele, 0);
            state->raw_matrix[index + 1] = founderallele_run(state, locus, personid, allele, 1);
        }
        else {
        */
            tmp = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, locus, allele));
            
            tmp2 = DESCENTGRAPH_GET(dg, 
                    DESCENTGRAPH_OFFSET(dg, last_personid, locus, last_allele));
            
            tmp3 = 1 - tmp; // XXX if i just do this inline i get a invalid write of 8 bytes!
            state->raw_matrix[index + tmp] = state->raw_matrix[index + tmp2];
            state->raw_matrix[index + tmp3] = founderallele_run(state, locus, personid, allele, tmp3);
        /*
        }
        */
    }
}

__global__ void msampler_sampling_kernel(struct gpu_state* state, int meiosis) {
    int personid = (state->founderallele_count / 2) + (meiosis / 2);
    int allele = meiosis % 2;
    struct geneticmap* map = GET_MAP(state);
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    int i, j;
    float total;
    
    __shared__ float sh_theta[1024];
    __shared__ float sh_inversetheta[1024];
    __shared__ float sh_matrix[1024][2];
    __shared__ int sh_descentgraph[1024][2];
    
    
    map_length = map->map_length;
    
    // we just have one block for now, 512 threads
    for(i = threadIdx.x; i < map_length; i += 256) {
        sh_descentgraph[i][0] = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, i, GPU_MATERNAL_ALLELE));
        sh_descentgraph[i][1] = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, i, GPU_PATERNAL_ALLELE));
        
        sh_matrix[i][0] = state->raw_matrix[i * 2];
        sh_matrix[i][1] = state->raw_matrix[(i * 2) + 1];
        
        if(i < (map_length - 1)) {
            sh_theta[i] = (float) MAP_THETA(map, i);
            sh_inversetheta[i] = (float) MAP_INVERSETHETA(map, i);
        }
    }
    
    __syncthreads();
    
    
    if(threadIdx.x == 0) {
        total = sh_matrix[0][0] + sh_matrix[0][1];
        sh_matrix[0][0] /= total;
        sh_matrix[0][1] /= total;
    
        // forward
        for(i = 1; i < map_length; ++i) {
            for(j = 0; j < 2; ++j) {
                sh_matrix[i][j] *= ( \
                    (sh_matrix[i-1][j]   * sh_inversetheta[i-1]) + \
                    (sh_matrix[i-1][1-j] * sh_theta[i-1]) \
                  );
            }
            
            total = sh_matrix[i][0] + sh_matrix[i][1];
            sh_matrix[i][0] /= total;
            sh_matrix[i][1] /= total;
        }
        
        // backward
        i = map_length - 1;
        sh_descentgraph[i][allele] = founderallele_sample2(state, sh_matrix[i][0], sh_matrix[i][1]);
        
        while(--i >= 0) {
            for(j = 0; j < 2; ++j) {
                sh_matrix[i][j] *= ((sh_descentgraph[i+1][allele] != j) ? sh_theta[i] : sh_inversetheta[i]);
            }
            
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
__global__ void msampler_window_sampling_kernel(struct gpu_state* state, int meiosis, int window_length) {
    int starting_locus = blockIdx.x * window_length;
    int current_locus;
    int personid = (state->founderallele_count / 2) + (meiosis / 2);
    int allele = meiosis % 2;
    struct geneticmap* map = GET_MAP(state);
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    int i, j;
    float total;
    int last_index;
    
    // these are all 20 because shared mem needs to be 
    // declared with a constant
    __shared__ float sh_theta[32];
    __shared__ float sh_inversetheta[32];
    __shared__ float sh_matrix[32][2];
    __shared__ int sh_descentgraph[32][2];
    
    
    map_length = map->map_length;
    
    current_locus = starting_locus + threadIdx.x;
    
    if((map_length - current_locus) < 2)
        return;
    
    if((threadIdx.x < window_length) && (current_locus < map_length)) {
        
        sh_descentgraph[threadIdx.x][0] = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, current_locus, GPU_MATERNAL_ALLELE));
        sh_descentgraph[threadIdx.x][1] = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, current_locus, GPU_PATERNAL_ALLELE));
        
        sh_matrix[threadIdx.x][0] = state->raw_matrix[current_locus * 2];
        sh_matrix[threadIdx.x][1] = state->raw_matrix[(current_locus * 2) + 1];
        
        if(current_locus < (map_length - 1)) {
            sh_theta[threadIdx.x]        = (float) MAP_THETA(       map, current_locus);
            sh_inversetheta[threadIdx.x] = (float) MAP_INVERSETHETA(map, current_locus);
        }
    }
    
    __syncthreads();
    
    
    if(threadIdx.x == 0) {
        total = sh_matrix[0][0] + sh_matrix[0][1];
        sh_matrix[0][0] /= total;
        sh_matrix[0][1] /= total;
        
        // forward
        for(i = 1; i < window_length; ++i) {
            if((starting_locus + i) >= map_length)
                break;
            
            for(j = 0; j < 2; ++j) {
                sh_matrix[i][j] *= ( \
                    (sh_matrix[i-1][j]   * sh_inversetheta[i-1]) + \
                    (sh_matrix[i-1][1-j] * sh_theta[i-1]) \
                  );
            }
            
            total = sh_matrix[i][0] + sh_matrix[i][1];
            sh_matrix[i][0] /= total;
            sh_matrix[i][1] /= total;
        
            last_index = i;
        }
        
        
        
        // backward
        i = last_index;
        sh_descentgraph[i][allele] = founderallele_sample2(state, sh_matrix[i][0], sh_matrix[i][1]);
        
        while(--i >= 0) {
            for(j = 0; j < 2; ++j) {
                sh_matrix[i][j] *= ((sh_descentgraph[i+1][allele] != j) ? sh_theta[i] : sh_inversetheta[i]);
            }
            
            sh_descentgraph[i][allele] = founderallele_sample2(state, sh_matrix[i][0], sh_matrix[i][1]);
        }
    }
    
    __syncthreads();
    
    
    if((threadIdx.x < window_length) && (current_locus < map_length)) {
        DESCENTGRAPH_SET(dg, DESCENTGRAPH_OFFSET(dg, personid, current_locus, GPU_MATERNAL_ALLELE), sh_descentgraph[threadIdx.x][0]);
        DESCENTGRAPH_SET(dg, DESCENTGRAPH_OFFSET(dg, personid, current_locus, GPU_PATERNAL_ALLELE), sh_descentgraph[threadIdx.x][1]);
    }
}
*/

__global__ void msampler_window_sampling_kernel(struct gpu_state* state, int meiosis, int window_length) {
    int starting_locus = blockIdx.x * window_length;
    int current_locus;
    int personid = (state->founderallele_count / 2) + (meiosis / 2);
    int allele = meiosis % 2;
    struct geneticmap* map = GET_MAP(state);
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    int i, j;
    float total;
    int last_index;
    
    // these are all 20 because shared mem needs to be 
    // declared with a constant
    __shared__ float sh_theta[32];
    __shared__ float sh_inversetheta[32];
    __shared__ float sh_matrix[32][2];
    __shared__ int sh_descentgraph[32][2];
    
    int buffer = 6;
    int left_buffer = starting_locus == 0 ? 0 : buffer;
    
    int window_starting_locus = starting_locus - left_buffer;
    int buffer_window_length = window_length + left_buffer + buffer;
    
    map_length = map->map_length;
    
    current_locus = window_starting_locus + threadIdx.x;
    
    if((map_length - current_locus) < 2)
        return;
    
    if((threadIdx.x < buffer_window_length) && (current_locus < map_length)) {
        
        sh_descentgraph[threadIdx.x][0] = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, current_locus, GPU_MATERNAL_ALLELE));
        sh_descentgraph[threadIdx.x][1] = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, current_locus, GPU_PATERNAL_ALLELE));
        
        sh_matrix[threadIdx.x][0] = state->raw_matrix[current_locus * 2];
        sh_matrix[threadIdx.x][1] = state->raw_matrix[(current_locus * 2) + 1];
        
        if(current_locus < (map_length - 1)) {
            sh_theta[threadIdx.x]        = (float) MAP_THETA(       map, current_locus);
            sh_inversetheta[threadIdx.x] = (float) MAP_INVERSETHETA(map, current_locus);
        }
    }
    
    __syncthreads();
    
    
    if(threadIdx.x == 0) {
        total = sh_matrix[0][0] + sh_matrix[0][1];
        sh_matrix[0][0] /= total;
        sh_matrix[0][1] /= total;
        
        // forward
        for(i = 1; i < buffer_window_length; ++i) {
            if((window_starting_locus + i) >= map_length)
                break;
            
            for(j = 0; j < 2; ++j) {
                sh_matrix[i][j] *= ( \
                    (sh_matrix[i-1][j]   * sh_inversetheta[i-1]) + \
                    (sh_matrix[i-1][1-j] * sh_theta[i-1]) \
                  );
            }
            
            total = sh_matrix[i][0] + sh_matrix[i][1];
            sh_matrix[i][0] /= total;
            sh_matrix[i][1] /= total;
        
            last_index = i;
        }
        
        
        
        // backward
        i = last_index;
        sh_descentgraph[i][allele] = founderallele_sample2(state, sh_matrix[i][0], sh_matrix[i][1]);
        
        while(--i >= 0) {
            for(j = 0; j < 2; ++j) {
                sh_matrix[i][j] *= ((sh_descentgraph[i+1][allele] != j) ? sh_theta[i] : sh_inversetheta[i]);
            }
            
            sh_descentgraph[i][allele] = founderallele_sample2(state, sh_matrix[i][0], sh_matrix[i][1]);
        }
    }
    
    __syncthreads();
    
    
    if((threadIdx.x >= left_buffer) && (threadIdx.x < (window_length + left_buffer)) && (current_locus < map_length)) {
        DESCENTGRAPH_SET(dg, DESCENTGRAPH_OFFSET(dg, personid, current_locus, GPU_MATERNAL_ALLELE), sh_descentgraph[threadIdx.x][0]);
        DESCENTGRAPH_SET(dg, DESCENTGRAPH_OFFSET(dg, personid, current_locus, GPU_PATERNAL_ALLELE), sh_descentgraph[threadIdx.x][1]);
    }
}

void run_gpu_msampler_reset_kernel(int numblocks, int numthreads, struct gpu_state* state, int meiosis, size_t shared) {
    msampler_reset_kernel<<<numblocks, numthreads, shared>>>(state, meiosis);
}

void run_gpu_msampler_likelihood_kernel(int numblocks, int numthreads, struct gpu_state* state, int meiosis, int last_meiosis, size_t shared) {
    msampler_likelihood_kernel<<<numblocks, numthreads, shared>>>(state, meiosis, last_meiosis);
}

void run_gpu_msampler_sampling_kernel(struct gpu_state* state, int meiosis) {
    msampler_sampling_kernel<<<1, 256>>>(state, meiosis);
}

void run_gpu_msampler_window_sampling_kernel(int numblocks, int numthreads, struct gpu_state* state, int meiosis, int window_length) {
    msampler_window_sampling_kernel<<<numblocks, numthreads>>>(state, meiosis, window_length);
}

void setup_msampler_kernel() {
    cudaFuncSetCacheConfig(msampler_likelihood_kernel, cudaFuncCachePreferShared);
    cudaFuncSetCacheConfig(msampler_sampling_kernel,   cudaFuncCachePreferShared);
}

