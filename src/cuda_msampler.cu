#include <stdio.h>
#include <float.h>

#include "cuda_common.h"


double TMP_LOG_ZERO = -DBL_MAX;

void founderallelegraph_print(struct gpu_state* state, int locus) {
    struct founderallelegraph* fag = GET_FOUNDERALLELEGRAPH(state, locus);
    struct adjacent_node* tmp;
    int i, j;
    
    printf("\nFOUNDER ALLELE GRAPH:\n");
    
	for(i = 0; i < state->founderallele_count; ++i) {
        printf("%d: ", i);
        for(j = 0; j < fag->num_neighbours[i]; ++j) {
            tmp = &(fag->graph[i][j]);
            printf("{%d, %d} ", tmp->id, tmp->label);
        }
        printf("\n");
    }
    printf("\n");
}

void print_descentgraph2(struct descentgraph* dg, int ped_length, int map_length) {
    int i, j;
    
    for(i = 0; i < ped_length; ++i) {
        printf("\t%d:\t", i);
        for(j = 0; j < map_length; ++j) {
            printf( "%d%d ",
                    DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, i, j, GPU_MATERNAL_ALLELE)),
                    DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, i, j, GPU_PATERNAL_ALLELE))
            );
        }
        printf("\n");
    }
    printf("\n");
}

double tmp_log_product(double a, double b) {
    return ((a == TMP_LOG_ZERO) or (b == TMP_LOG_ZERO)) ? TMP_LOG_ZERO : a + b;
}

double tmp_log_sum(double a, double b) {
    if(a == TMP_LOG_ZERO)
        return b;
    
    if(b == TMP_LOG_ZERO)
        return a;
    
    return log(exp(b - a) + 1) + a;
}

int founderallele_add(struct founderallelegraph* fag, int mat_fa, int pat_fa, int g) {
    struct adjacent_node* tmp;
    int i;
    int neighbours;
    
    if(g == GPU_GENOTYPE_UNTYPED) {
        return 1;
    }
    
    // check to see if it exists
    // if it does, then the edge needs to have the same label
    neighbours = fag->num_neighbours[mat_fa];
    
    for(i = 0; i < neighbours; ++i) {
        tmp = &(fag->graph[mat_fa][i]);
        
        if(tmp->id == pat_fa) {
            return tmp->label == g;
        }
    }
    
    tmp = &(fag->graph[mat_fa][fag->num_neighbours[mat_fa]]);
    tmp->id = pat_fa;
    tmp->label = g;
    fag->num_neighbours[mat_fa]++;
    
    if(mat_fa != pat_fa) {
        tmp = &(fag->graph[pat_fa][fag->num_neighbours[pat_fa]]);
        tmp->id = mat_fa;
        tmp->label = g;
        fag->num_neighbours[pat_fa]++;
    }
        
    return 1;
}

int founderallelegraph_populate(struct gpu_state* state, int locus) {
	struct person* p;
	struct founderallelegraph* fag = GET_FOUNDERALLELEGRAPH(state, locus);
	struct descentgraph* dg = GET_DESCENTGRAPH(state);
	int g;
	int i;
	int legal = 1;
	int pid;
	int founderalleles[256];
	int parent_allele;
	int mat, pat;
	
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
	        
	        if(! founderallele_add(fag, mat, pat, g)) {
                legal = 0;
                printf("illegal (locus = %d, person = %d, [%d %d %d])!\n", locus, i, mat, pat, g);
            }
        }
	}
	
	return legal;
}

// for loops in the allele graph, hetero is always a contradiction
int correct_alleles_loop(int g, int allele1) {
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

int correct_alleles(int g, int allele1) {
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

int correct_alleles(int g, int allele1, int allele2) {
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


int legal(struct gpu_state* state, int locus, int* component, int clength, int* q, int qlength) {
//bool FounderAlleleGraph2::legal(GraphComponent& gc, vector<unsigned>& assignment) {
    //int node = gc[assignment.size() - 1];
    //int allele = assignment.back();
    int node = component[qlength-1];
    int allele = q[qlength-1];
    int i, j;  
    
    //AdjacencyRecord& tmp = matrix[node];
    struct founderallelegraph* fag = GET_FOUNDERALLELEGRAPH(state, locus);
    struct adjacent_node* adj;

    //for(unsigned i = 0; i < tmp.size(); ++i) {
    for(i = 0; i < fag->num_neighbours[node]; ++i) {
        //FounderAlleleNode& adj = tmp[i];
        adj = &(fag->graph[node][i]);
        
        // if there is a loop
        //if(adj.id == node) {
        //    if(not correct_alleles_loop(adj.label, allele)) {
        //        return false;
        //    }
        //    continue;
        //}
        
        // if there is a loop
        if(adj->id == node) {
            if(! correct_alleles_loop(adj->label, allele)) {
                return 0;
            }
            continue;
        }
        
        // test if compatible with label
        //if(not correct_alleles(adj.label, allele)) {
        //    return false;
        //}
        
        // test if compatible with label
        if(! correct_alleles(adj->label, allele)) {
            return 0;
        }
        
        // find offset of adjacent node in assignment vector
        //unsigned j;
        //for(j = 0; j < gc.size(); ++j) {
        //    if(gc[j] == adj.id)
        //        break;
        //}
        
        // find offset of adjacent node in assignment vector
        for(j = 0; j < clength; ++j) {
            if(component[j] == adj->id) {
                break;
            }
        }
        
        // error if not found
        //if(j == gc.size()) {
        //    fprintf(stderr, "Error: an adjacent allele in the graph was not "
        //                    "found in the same component (%s:%d)", 
        //                    __FILE__, __LINE__);
        //    abort();
        //}

        // error if not found
        if(j == clength) {
            // XXX abort on gpu?
            abort();
            //return 0;
        }
        
        // if not assigned yet, then ignore
        //if(j > (assignment.size() - 1))
        //    continue;

        // if not assigned yet, then ignore
        if(j > (qlength-1))
            continue;
        
        // if assigned, then test legality
        //if(not correct_alleles(adj.label, allele, assignment[j])) {
        //    return false;
        //}

        // if assigned, then test legality
        if(! correct_alleles(adj->label, allele, q[j])) {
            return 0;
        }

    }
    
    return 1;
}

double component_likelihood(struct gpu_state* state, int locus, int* q, int length) {
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

double founderallelegraph_enumerate(struct gpu_state* state, int locus, int* component, int cindex) {
    int q[128];
    int qindex = 0;
    int skip = 0;
    double prob = 0.0;
    int i;

    // <debug>
    printf("* ");
    for(i = 0; i < cindex; ++i)
        printf("%d ", component[i]);
    printf("\n");
    // </debug>
    
    if(cindex == 1) {
        return 1.0;
    }

    while(1) {
        
        while(qindex != cindex) {
            q[qindex++] = 1; // push
            if(! legal(state, locus, component, cindex, q, qindex)) {
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
            
            if(legal(state, locus, component, cindex, q, qindex)) {
                skip = 0;
                break;
            }
        }
    }
       
no_more_assignments:
    
    return prob;
}

double founderallelegraph_likelihood(struct gpu_state* state, int locus) {
    int q[128];
    int qindex = 0;
    
    int component[128];
    int cindex;
    
    int visited[128];
    int total = state->founderallele_count;
    
    int i;
    int tmp;
    
    double tmp_prob;
    double prob = 0.0;
    
    struct founderallelegraph* fag = GET_FOUNDERALLELEGRAPH(state, locus);
    struct adjacent_node* tmp2;
    
    
    for(i = 0; i < 128; ++i) {
        visited[i] = _WHITE;
    }
    
    do {
        cindex = 0;
        
        // find start point
        for(i = 0; i < state->founderallele_count; ++i) {
            if(visited[i] == _WHITE) {
                visited[i] = _GREY;
                q[qindex++] = i;
                break;
            }
        }
        
        while(qindex != 0) {
            tmp = q[--qindex];
            
            for(i = 0; i < fag->num_neighbours[tmp]; ++i) {
                tmp2 = &(fag->graph[tmp][i]);
                
                if(visited[tmp2->id] == _WHITE) {
                    visited[tmp2->id] = _GREY;
                    q[qindex++] = tmp2->id;
                }
            }
            
            visited[tmp] = _BLACK;
            total--;
            
            component[cindex++] = tmp;
        }
        
        tmp_prob = founderallelegraph_enumerate(state, locus, component, cindex);
		tmp_prob = ((tmp_prob == 0.0) ? TMP_LOG_ZERO : log(tmp_prob));
		
		prob = tmp_log_product(tmp_prob, prob);
    
    } while(total != 0);
    
    return prob;
}

int founderallele_sample(struct founderallelegraph* fag) {
    double total; // = fag->prob[0] + fag->prob[1];
    
    if((fag->prob[0] == TMP_LOG_ZERO) && (fag->prob[1] == TMP_LOG_ZERO)) {
        abort();
    }
    
    if(fag->prob[0] == TMP_LOG_ZERO)
        return 1;
        
    if(fag->prob[1] == TMP_LOG_ZERO)
        return 0;
    
    total = fag->prob[0] + fag->prob[1];
    
    return (rand() / double(RAND_MAX)) < (fag->prob[0] / total) ? 0 : 1;
}

double founderallele_run(struct gpu_state* state, int locus, int personid, int allele, int value) {
    struct founderallelegraph* fag = GET_FOUNDERALLELEGRAPH(state, locus);
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    int i;
    int tmp;
    int populate_legal;
    double prob = TMP_LOG_ZERO;
    
    // save the previous value
    i = DESCENTGRAPH_OFFSET(dg, personid, locus, allele);
    tmp = DESCENTGRAPH_GET(dg, i);
    DESCENTGRAPH_SET(dg, i, value);
    
    for(i = 0; i < state->founderallele_count; ++i) {
        fag->num_neighbours[i] = 0;
    }
    
    populate_legal = founderallelegraph_populate(state, locus);
    
    if(populate_legal) {
        prob = founderallelegraph_likelihood(state, locus);
    }
    
    //printf("%f\n", prob);
    
    // restore previous value
    i = DESCENTGRAPH_OFFSET(dg, personid, locus, allele);
    DESCENTGRAPH_SET(dg, i, tmp);
    
    return prob;
}

void msampler_kernel(struct gpu_state* state, int meiosis) {
    //int locus = blockIdx.x;
    //sampler_run(state, locus);
    
    /*
    printf("fa_sequence = ");
    for(i = 0; i < state->pedigree_length; ++i) {
        printf("%d ", state->fa_sequence[i]);
    }
    printf("\n");
    */
    
    int locus;
    int personid = (state->founderallele_count / 2) + (meiosis / 2);
    int allele = meiosis % 2;
    struct founderallelegraph* fag;
    struct geneticmap* map;
    struct descentgraph* dg = GET_DESCENTGRAPH(state);
    int i, j;
    
    map = GET_MAP(state);
        
    // build fags
    for(locus = 0; locus < map->map_length; ++locus) {
        fag = GET_FOUNDERALLELEGRAPH(state, locus);
        
        fag->prob[0] = founderallele_run(state, locus, personid, allele, 0);
        fag->prob[1] = founderallele_run(state, locus, personid, allele, 1);
    }
    
    
    // sampling
    
    // forward
    for(i = 1; i < state->map->map_length; ++i) {
        for(j = 0; j < 2; ++j) {
            state->graphs[i].prob[j] = tmp_log_product( \
                                            state->graphs[i].prob[j], \
                                            tmp_log_sum( \
                                                tmp_log_product(state->graphs[i-1].prob[j],   MAP_THETA(map, i-1)), \
                                                tmp_log_product(state->graphs[i-1].prob[1-j], MAP_INVERSETHETA(map, i-1)) \
                                            ) \
                                        );
        }
    }
    
    /*
    // print
    for(locus = 0; locus < state->map->map_length; ++locus) {
        fag = GET_FOUNDERALLELEGRAPH(state, locus);
        printf("%d %.3f %.3f\n", locus, fag->prob[0], fag->prob[1]);
    }
    */
    
    // backward
    i = map->map_length - 1;
    fag = GET_FOUNDERALLELEGRAPH(state, i);
    DESCENTGRAPH_SET(dg, DESCENTGRAPH_OFFSET(dg, personid, i, allele), founderallele_sample(fag));
    
    while(--i >= 0) {
        fag = GET_FOUNDERALLELEGRAPH(state, i);
        
        for(j = 0; j < 2; ++j) {
            fag->prob[j] = tmp_log_product(fag->prob[j], ((DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, personid, i+1, allele)) != j) ? \
                log(MAP_THETA(map, i)) : \
                log(MAP_INVERSETHETA(map, i))));
        }
        
        DESCENTGRAPH_SET(dg, DESCENTGRAPH_OFFSET(dg, personid, i, allele), founderallele_sample(fag));
    }
    
    /*
    // print
    for(locus = 0; locus < state->map->map_length; ++locus) {
        fag = GET_FOUNDERALLELEGRAPH(state, locus);
        printf("%d %.3f %.3f\n", locus, fag->prob[0], fag->prob[1]);
    }
    */
}

void run_gpu_msampler_kernel(int numblocks, int numthreads, struct gpu_state* state, int meiosis) {
    //msampler_kernel<<<numblocks, numthreads>>>(state, meiosis);
    msampler_kernel(state, meiosis);
}

