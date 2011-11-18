#include "cuda_common.h"

__device__ void print_ints(int* data, int length) {
    int i;
    
    for(i = 0; i < length; ++i) {
        printf("%d ", data[i]);
    }
    
    printf("\n");
}

__device__ void print_rfunction(struct rfunction* r) {
    printf( "\tid: %d\n"
            "\ttype: %d\n"
            "\tprev1: %d\n"
            "\tprev2: %d\n",
            r->id,
            r->peel_type,
            r->prev1 != NULL ? r->prev1->id : -1,
            r->prev2 != NULL ? r->prev2->id : -1
    );
    
    printf("\tcutset: ");
    print_ints(r->cutset, r->cutset_length);
    printf("\n");
}

__device__ void print_person(struct person* p) {
    printf( "\tid: %d\n"
            "\tmother: %d\n"
            "\tfather: %d\n"
            "\tfounder: %d\n"
            "\ttyped: %d\n"
            "\tprobs: %.2e %.2e %.2e %.2e\n",
            p->id,
            p->mother,
            p->father,
            p->isfounder,
            p->istyped,
            p->prob[0],
            p->prob[1],
            p->prob[2],
            p->prob[3]
    );
    
    printf("\tgenotypes: ");
    print_ints(p->genotypes, p->genotypes_length);
    printf("\n");
}

__device__ void print_map(struct geneticmap* map) {
    int i, j;
    
    for(i = 0; i < (map->map_length - 1); ++i) {
        printf("\t%d\t%.3f %.3f\n", i, map->thetas[i], map->inversethetas[i]);
    }
    printf("\n");
    
    for(i = 0; i < map->map_length; ++i) {
        printf("\t%d:\t", i);
        for(j = 0; j < 4; ++j) {
            printf("%.3f ", MAP_PROB(map, i, j));
        }
        printf("\n");
    }
    printf("\n");
}

__device__ void print_descentgraph(struct descentgraph* dg, int ped_length, int map_length) {
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

__device__ void print_everything(struct gpu_state* state) {
    int i;
    
    printf("RFUNCTIONS:\n");
    
    for(i = 0; i < state->functions_length; ++i) {
        print_rfunction(&state->functions[i]);
    }
    printf("\n");
    
    printf("PEDIGREE:\n");
    for(i = 0; i < state->pedigree_length; ++i) {
        print_person(&state->pedigree[i]);
    }
    printf("\n");
    
    printf("MAP:\n");
    print_map(state->map);
    printf("\n");
    
    printf("DESCENTGRAPH:\n");
    print_descentgraph(state->dg, state->pedigree_length, state->map->map_length);
    printf("\n");
}

__global__ void print_pedigree(struct person* p, int length) {
    int i;
    
    printf("device: sizeof(struct person) = %d\n", sizeof(struct person));
    
    for(i = 0; i < length; ++i) {
        print_person(&p[i]);
    }
}

__global__ void print_kernel(struct gpu_state* state) {
    if(threadIdx.x == 0)
        print_everything(state);
}

void run_gpu_print_pedigree_kernel(struct person* p, int length) {
    print_pedigree<<<1, 1>>>(p, length);
}

void run_gpu_print_kernel(struct gpu_state* state) {
    print_kernel<<<1, 1>>>(state);
}

