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

int get_founderallele(struct gpu_state* state, int person, int locus, int allele) {
    int current = person;
    int parent_allele = allele;
	struct person* p;
	struct descentgraph* dg = GET_DESCENTGRAPH(state);
    
    while(1) {
		p = GET_PERSON(state, current);

        if(PERSON_ISFOUNDER(p)) {
            return (current * 2) + parent_allele;
        }
        
        if(parent_allele == GPU_PATERNAL_ALLELE) {
            parent_allele = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, current, locus, parent_allele));
            current = PERSON_FATHER(p);
        }
        else {
            parent_allele = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(dg, current, locus, parent_allele));
            current = PERSON_MOTHER(p);
		}
    }
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
	int mat;
	int pat;
	int g;
	int i;
	int legal = 1;
    
	for(i = 0; i < state->pedigree_length; ++i) {
        p = GET_PERSON(state, i);
        
        if(PERSON_ISTYPED(p)) {
            mat = get_founderallele(state, i, locus, GPU_MATERNAL_ALLELE);
            pat = get_founderallele(state, i, locus, GPU_PATERNAL_ALLELE);
            g = PERSON_GENOTYPE(p, locus);
            
            if(!founderallele_add(fag, mat, pat, g)) {
                legal = 0;
                printf("illegal (locus = %d, person = %d, [%d %d %d])!\n", locus, i, mat, pat, g);
            }
		}
        
		//__syncthreads();
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

void msampler_run(struct gpu_state* state, int locus) {
    struct founderallelegraph* fag = GET_FOUNDERALLELEGRAPH(state, locus);
    int i;
    double prob = TMP_LOG_ZERO;
    
    for(i = 0; i < state->founderallele_count; ++i) {
        fag->num_neighbours[i] = 0;
    }
    
    if(founderallelegraph_populate(state, locus)) {
        prob = founderallelegraph_likelihood(state, locus);
    }
    
    printf("%f\n", prob);
    
    founderallelegraph_print(state, locus);
    print_descentgraph2(GET_DESCENTGRAPH(state), state->pedigree_length, state->map->map_length);
    
    //__syncthreads();
    
    // enumerate + likelihood
    
}

void msampler_kernel(struct gpu_state* state) {
    //int locus = blockIdx.x;
    //sampler_run(state, locus);
    
    int locus;
    
    // build fags
    for(locus = 0; locus < state->map->map_length; ++locus) {
        msampler_run(state, locus);
    }
    
    printf("\n\n");
    
    abort();
    
    // sample
    // XXX
}

void run_gpu_msampler_kernel(int numblocks, int numthreads, struct gpu_state* state) {
    //msampler_kernel<<<numblocks, numthreads>>>(state);
    msampler_kernel(state);
}

