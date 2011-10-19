#include "cuda_common.h"


void founderallelegraph_print(struct gpu_state* state, int locus) {
    struct founder_allele_graph* fag = GET_FOUNDERALLELEGRAPH(state, locus);
    struct adjacency_node* tmp;
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
        
        parent_allele = DESCENTGRAPH_GET(dg, DESCENTGRAPH_OFFSET(current, locus, parent_allele));
                
        if(parent_allele == GPU_PATERNAL_ALLELE) {
            current = PERSON_FATHER(p);
        }
        else {
            current = PERSON_MOTHER(p);
		}
    }
}

int founderallele_add(struct founder_allele_graph* fag, int mat_fa, int pat_fa, int g) {
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
	struct founder_allele_graph* fag = GET_FOUNDERALLELEGRAPH(state, locus);
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
            g = PERSON_GENOTYPE(locus);
            
            if(!founderallelegraph_add(fag, mat, pat, g)) {
                legal = 0;
            }
		}
        
		//__syncthreads();
	}
	
	return legal;
}

void msampler_run(struct gpu_state* state, int locus) {
    int i;
    struct founder_allele_graph* graph = GET_FOUNDERALLELEGRAPH(state, locus);
    struct adjacency_node* nodes;
    int num_neighbours = 0;
    
    // construct
    founderallelegraph_populate(state, locus);
    founderallelegraph_print(state, locus);
    
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
    
    // sample
    // XXX
}

void run_gpu_msampler_kernel(int numblocks, int numthreads, struct gpu_state* state) {
    //msampler_kernel<<<numblocks, numthreads>>>(state);
    msampler_kernel(state);
}

