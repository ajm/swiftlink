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
	
    //int tmp = (locus % 8) * 16 * state->founderallele_count;
    int tmp = (14 * state->founderallele_count) + (2 * state->pedigree_length);
    tmp += (tmp % 4);
    tmp *= (locus % 8);
    
    /*
     int8_t* group_membership   = ( int8_t*) &extern_pool[tmp];
     int8_t* group_fixed        = ( int8_t*) &extern_pool[tmp +      state->founderallele_count];
    uint8_t* group_active       = (uint8_t*) &extern_pool[tmp + (2 * state->founderallele_count)];
    uint8_t* group_size         = (uint8_t*) &extern_pool[tmp + (3 * state->founderallele_count)];
    uint8_t* allele_assignment  = (uint8_t*) &extern_pool[tmp + (4 * state->founderallele_count)];
    float*   prob               = ( float* ) &extern_pool[tmp + (6 * state->founderallele_count)];
    uint8_t* edge_list          = (uint8_t*) &extern_pool[tmp + (14 * state->founderallele_count)];
    */
    float*   prob               = ( float* ) &extern_pool[tmp];
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
    
    
    //printf("\nlocus = %d\n\n", locus);
    //printf("\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n\n", group_membership, group_fixed, group_active, group_size, edge_list, allele_assignment, prob);
    
    
    for(i = 0; i < state->founderallele_count; ++i) {
        group_membership[i] = GPU_DEFAULT_COMPONENT;
        group_fixed[i] = -1;
        group_active[i] = 0;
        //allele_assignment[ind(0, i)] = 9;
        //allele_assignment[ind(1, i)] = 9;
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
	    
	    //printf("* %d (%d - %d) %d\n", pid, edge_list[pid * 2], edge_list[(pid * 2) + 1], PERSON_GENOTYPE(p, locus));
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
        
        tmp = i * 2;
        mat_fa = edge_list[tmp];
        pat_fa = edge_list[++tmp];
        
        /*
        printf("group_membership: ");
        for(j = 0; j < state->founderallele_count; ++j) {
            printf("%d ", group_membership[j]);
        }
        printf("\n");
        
        printf("group_active: ");
        for(j = 0; j < state->founderallele_count; ++j) {
            printf("%d ", group_active[j]);
        }
        printf("\n");
        
        printf("group_fixed: ");
        for(j = 0; j < state->founderallele_count; ++j) {
            printf("%d ", group_fixed[j]);
        }
        printf("\n");
        
        printf("allele_assignment: ");
        for(j = 0; j < state->founderallele_count; ++j) {
            printf("%d ", allele_assignment[ind(0,j)]);
        }
        printf("\n");
        
        printf("allele_assignment: ");
        for(j = 0; j < state->founderallele_count; ++j) {
            printf("%d ", allele_assignment[ind(1,j)]);
        }
        printf("\n\n");
        */
        
        // autozygous
        if(mat_fa == pat_fa) {
            if(g == GPU_GENOTYPE_AB) {
                //printf("zp %s:%d\n", __FILE__, __LINE__);
                return 0.0;
            }
            
            group1 = group_membership[mat_fa];
            
            if(group1 != GPU_DEFAULT_COMPONENT) {
                // already belongs to a component
                fixed1 = group_fixed[group1];
                if(fixed1 != -1) {
                    if(g != allele_assignment[ind(fixed1, mat_fa)]) {
                        //printf("zp %s:%d\n", __FILE__, __LINE__);
                        return 0.0;
                    }
                }
                else {
                    if(allele_assignment[ind(0, mat_fa)] == g) {
                        group_fixed[group1] = 0;
                    }
                    else if(allele_assignment[ind(1, mat_fa)] == g) {
                        group_fixed[group1] = 1;
                    }
                    else {
                        //printf("zp %s:%d\n", __FILE__, __LINE__);
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
                prob[ind(0, group_index)] = (g == GPU_GENOTYPE_AA ? major_freq : minor_freq);
                
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
                                //printf("zp %s:%d\n", __FILE__, __LINE__);
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
                                }
                            }
                            else {
                                if(legal1) {
                                    group_fixed[group1] = 1;
                                }
                                else {
                                    /*
                                    printf("zp %s:%d %d [%d:%d:%d:%d:%d]\n", \
                                        __FILE__, __LINE__, i, g, \
                                        allele_assignment[ind(0, mat_fa)], \
                                        allele_assignment[ind(0, pat_fa)], \
                                        allele_assignment[ind(1, mat_fa)], \
                                        allele_assignment[ind(1, pat_fa)]);
                                    */
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
                                    //printf("zp %s:%d\n", __FILE__, __LINE__);
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
                                    //printf("zp %s:%d\n", __FILE__, __LINE__);
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
                                //printf("zp %s:%d\n", __FILE__, __LINE__);
                                return 0.0;
                            }
                        }
                        else {
                            // neither group1 nor group2 are fixed
                            legal0 = legal(g, allele_assignment[ind(0, mat_fa)], allele_assignment[ind(0, pat_fa)]);
                            legal1 = legal(g, allele_assignment[ind(1, mat_fa)], allele_assignment[ind(0, pat_fa)]);
                            legal2 = legal(g, allele_assignment[ind(0, mat_fa)], allele_assignment[ind(1, pat_fa)]);
                            legal3 = legal(g, allele_assignment[ind(1, mat_fa)], allele_assignment[ind(1, pat_fa)]);
                            
                            /*
                            printf("legal0 = %d: %d, %d, %d\n", legal0, g, allele_assignment[ind(0, mat_fa)], allele_assignment[ind(0, pat_fa)]);
                            printf("legal1 = %d: %d, %d, %d\n", legal1, g, allele_assignment[ind(1, mat_fa)], allele_assignment[ind(0, pat_fa)]);
                            printf("legal2 = %d: %d, %d, %d\n", legal2, g, allele_assignment[ind(0, mat_fa)], allele_assignment[ind(1, pat_fa)]);
                            printf("legal3 = %d: %d, %d, %d\n", legal3, g, allele_assignment[ind(1, mat_fa)], allele_assignment[ind(1, pat_fa)]);
                            */
                            
                            if(! (legal0 || legal1 || legal2 || legal3)) {
                                //printf("zp %s:%d\n", __FILE__, __LINE__);
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
                        int tmp;
                        
                        if(fixed1 != fixed2) {
                            // XXX for now, always keep group1
                            group_fixed[group1] = fixed1;
                            //combine_components(group1, group2, 1);
                        }
                        else {
                            if(fixed1 == -2) {
                                group_fixed[group1] = -1;
                                //combine_components(group1, group2, 1);
                            }
                            else {
                                group_fixed[group1] = fixed1;
                                //combine_components(group1, group2, 0);
                                flip = 0;
                            }
                        }
                        
                        // begin combine component code
                        for(j = 0; j < state->founderallele_count; ++j) {
                            if(group_membership[j] == group2) {
                                group_membership[j] = group1;
                                
                                if(flip) {
                                    tmp = allele_assignment[ind(0, j)];
                                    allele_assignment[ind(0, j)] = allele_assignment[ind(1, j)];
                                    allele_assignment[ind(1, j)] = tmp;
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
                        
                        //printf("combine i=%d flip=%d (%d, %d, %d, %d)\n", i, flip, legal0, legal1, legal2, legal3);
                        //printf("group1 = %d\ngroup2 = %d\n", group1, group2);
                        
                        group_size[group1] += group_size[group2];
                        group_active[group2] = 0; // false
                        // end combine component code
                        
                        --num_groups;
                    }
                }
                else {
                    if(fixed1 != -1) {
                        int tmp = get_other_allele(g, allele_assignment[ind(fixed1, mat_fa)]);
                        
                        if(tmp == GPU_GENOTYPE_UNTYPED) {
                            //printf("zp %s:%d\n", __FILE__, __LINE__);
                            return 0.0;
                        }
                        else {
                            allele_assignment[ind(fixed1, pat_fa)] = tmp;
                            prob[ind(fixed1, group1)] *= (tmp == GPU_GENOTYPE_AA ? major_freq : minor_freq);
                        }
                    }
                    else {
                        int tmp0 = get_other_allele(g, allele_assignment[ind(0, mat_fa)]);
                        int tmp1 = get_other_allele(g, allele_assignment[ind(1, mat_fa)]);
                        
                        if(tmp0 != GPU_GENOTYPE_UNTYPED) {
                            if(tmp1 != GPU_GENOTYPE_UNTYPED) {
                                allele_assignment[ind(0, pat_fa)] = tmp0;
                                allele_assignment[ind(1, pat_fa)] = tmp1;
                                prob[ind(0, group1)] *= (tmp0 == GPU_GENOTYPE_AA ? major_freq : minor_freq);
                                prob[ind(1, group1)] *= (tmp1 == GPU_GENOTYPE_AA ? major_freq : minor_freq);
                            }
                            else {
                                allele_assignment[ind(0, pat_fa)] = tmp0;
                                prob[ind(0, group1)] *= (tmp0 == GPU_GENOTYPE_AA ? major_freq : minor_freq);
                                group_fixed[group1] = 0;
                            }
                        }
                        else {
                            if(tmp1 != GPU_GENOTYPE_UNTYPED) {
                                allele_assignment[ind(1, pat_fa)] = tmp1;
                                prob[ind(1, group1)] *= (tmp1 == GPU_GENOTYPE_AA ? major_freq : minor_freq);
                                group_fixed[group1] = 1;
                            }
                            else {
                                //printf("zp %s:%d\n", __FILE__, __LINE__);
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
                    int tmp = get_other_allele(g, allele_assignment[ind(fixed2, pat_fa)]);
                    
                    if(tmp == GPU_GENOTYPE_UNTYPED) {
                        //printf("zp %s:%d\n", __FILE__, __LINE__);
                        return 0.0;
                    }
                    else {
                        allele_assignment[ind(fixed2, mat_fa)] = tmp;
                        prob[ind(fixed2, group2)] *= (tmp == GPU_GENOTYPE_AA ? major_freq : minor_freq);
                    }
                }
                else {
                    int tmp0 = get_other_allele(g, allele_assignment[ind(0, pat_fa)]);
                    int tmp1 = get_other_allele(g, allele_assignment[ind(1, pat_fa)]);
                    
                    if(tmp0 != GPU_GENOTYPE_UNTYPED) {
                        if(tmp1 != GPU_GENOTYPE_UNTYPED) {
                            allele_assignment[ind(0, mat_fa)] = tmp0;
                            allele_assignment[ind(1, mat_fa)] = tmp1;
                            prob[ind(0, group2)] *= (tmp0 == GPU_GENOTYPE_AA ? major_freq : minor_freq);
                            prob[ind(1, group2)] *= (tmp1 == GPU_GENOTYPE_AA ? major_freq : minor_freq);
                        }
                        else {
                            allele_assignment[ind(0, mat_fa)] = tmp0;
                            prob[ind(0, group2)] *= (tmp0 == GPU_GENOTYPE_AA ? major_freq : minor_freq);
                            group_fixed[group2] = 0;
                        }
                    }
                    else {
                        if(tmp1 != GPU_GENOTYPE_UNTYPED) {
                            allele_assignment[ind(1, mat_fa)] = tmp1;
                            prob[ind(1, group2)] *= (tmp1 == GPU_GENOTYPE_AA ? major_freq : minor_freq);
                            group_fixed[group2] = 1;
                        }
                        else {
                            //printf("zp %s:%d\n", __FILE__, __LINE__);
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
                    prob[ind(0, group_index)] = pow(g == GPU_GENOTYPE_AA ? major_freq : minor_freq, 2);
                    group_fixed[group_index] = 0;
                }
                
                // neither in a group
                group_membership[mat_fa] = group_membership[pat_fa] = group_index;
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
        if(group_active[i] && (group_size[i] > 1)) {
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
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    int i, tmp;
    double prob = _LOG_ZERO;
    
    // save the previous value
    i = DESCENTGRAPH_OFFSET(dg, personid, locus, allele);
    tmp = DESCENTGRAPH_GET(dg, i);
    DESCENTGRAPH_SET(dg, i, value);
    
    prob = founderallelegraph_likelihood(state, locus);
    //printf("prob = %e\n", prob);
    prob = (prob == 0.0 ? _LOG_ZERO : logf(prob));
    
    // restore previous value
    i = DESCENTGRAPH_OFFSET(dg, personid, locus, allele);
    DESCENTGRAPH_SET(dg, i, tmp);
    
    return prob;
}

__global__ void msampler_likelihood_kernel(struct gpu_state* state, int meiosis) {
    //int locus = (blockIdx.x * NUM_THREADS) + threadIdx.x;
    int locus = ((blockIdx.x * 256) + threadIdx.x) / 32;
    int personid = (state->founderallele_count / 2) + (meiosis / 2);
    int allele = meiosis % 2;
    struct founderallelegraph* fag;
    struct geneticmap* map = GET_MAP(state);
    
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    int tmp, tmp2, tmp3;
    //int i ;
    
    //if(threadIdx.x == 0) {
    //for(i = 0; i < map->map_length; ++i) {
    //    locus = i;
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
    //}
    //}
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
