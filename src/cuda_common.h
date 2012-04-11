#ifndef GPU_RFUNCTION_H_
#define GPU_RFUNCTION_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cuda_quiet.h"

#include "tinymt/tinymt32_host.h"

//#define NUM_THREADS 96
#define NUM_ALLELES 4
#define GPU_DEFAULT_COMPONENT -1

enum {
    GPU_TRAIT_A,
    GPU_TRAIT_B
};

enum { 
    GPU_TRAIT_AA,
    GPU_TRAIT_BA,
    GPU_TRAIT_AB,
    GPU_TRAIT_BB
};

enum {
	GPU_GENOTYPE_UNTYPED,
	GPU_GENOTYPE_AB,
	GPU_GENOTYPE_AA,
	GPU_GENOTYPE_BB
};

enum {
    GPU_CHILD_PEEL,
    GPU_PARTNER_PEEL,
    GPU_PARENT_PEEL
};

enum {
    GPU_MATERNAL_ALLELE,
    GPU_PATERNAL_ALLELE
};

enum {
    _WHITE,
    _GREY,
    _BLACK
};

struct gpu_state {
    // lsampler stuff
    struct rfunction* functions;
    int functions_length;
    int functions_per_locus;
    
    // msampler stuff
    //struct founderallelegraph* graphs;
    double* raw_matrix;
    int founderallele_count;
    
    // constant information
    struct person* pedigree;
    int pedigree_length;
    
    struct geneticmap* map;
    struct descentgraph* dg;
    
    // scoring
    double* lodscores;
    
    // random number generation
    curandState* randstates;
    tinymt32_status_t* randmt;
    
    int* fa_sequence;
};

/*
struct adjacent_node {
    int id;
    int label;
};

struct adjacent_node2 {
    //short int id;
    //short int label;
    uint8_t id;
    uint8_t label;
};

struct founderallelegraph {
    float prob[2];
    //int* num_neighbours;            // length num founder alleles
    //struct adjacent_node* graph;   // num founder alleles * (num founder alleles + 1)
};

#define FAG_GET_NODE(fag_ptr, node)         (&(fag_ptr)->graph[(node)])
#define FAG_GET_NEIGHBOURS(fag_ptr, node)   (&(fag_ptr)->num_neighbours[(node)])
#define NODE_ID(node_ptr, n)                ((node_ptr)[(n)]->id)
#define NODE_LABEL(node_ptr, n)             ((node_ptr)[(n)]->label)
*/

struct geneticmap {
    double* thetas;
    double* inversethetas;
    double* halfthetas;
    double* halfinversethetas;
    double* markerprobs;
    double* allelefreqs;
    int map_length;
};

// THE 'peel_node' MUST ALWAYS BE LOCATED
// AT cutset[cutset_length - 1], SO IT CAN 
// BE IGNORED EASILY
struct rfunction {
    int id;

    int peel_type;
    int peel_node;
        
    int* cutset;                // eg: [1,2,3] (including the 'peel_node')
    int cutset_length;
    
    // how are these two matrices going to be treated?
    // could the index to one be a factor of the other
    // eg: presum_matrix[(x,y,z)] and matrix[(x,y,z) / z]?
    
    double* presum_matrix;
    int presum_length;          // 4 ** (cutset_length + 1)
    
    double* matrix;
    int matrix_length;          // 4 ** cutset_length
    
    double* transmission;       // 64 * children_length
    int* children;
    int children_length;
    
    struct rfunction* prev1;    // must be NULL if not used
    struct rfunction* prev2;    // must be NULL if not used
};

struct person {
    int id;
    int father;
    int mother;
    int isfounder;
    int istyped;
    
    int* genotypes;
    int genotypes_length;
            
    float prob[4];  // <--- XXX if this is double, then the sizeof differs between host and device XXX
                    //  this is probably due to the genotypes pointer which is different on different
                    //  architectures and maybe the compiler is trying to add some padding as things
                    //  might need to be word aligned in different situations
};

struct descentgraph {
    int* graph;
    int graph_length;
    int subgraph_length;
    int* padding;
    double transmission_prob;
};

#define GET_RFUNCTION(state_ptr, n, locus) (&(state_ptr)->functions[((locus) * (state_ptr)->functions_per_locus) + (n)])
#define GET_PERSON(state_ptr, n) (&(state_ptr)->pedigree[(n)])
#define GET_DESCENTGRAPH(state_ptr) ((state_ptr)->dg)
#define GET_MAP(state_ptr) ((state_ptr)->map)
#define GET_FOUNDERALLELEGRAPH(state_ptr, n) (&(state_ptr)->graphs[(n)])

#define RFUNCTION_GET(rf_ptr, index)        ((rf_ptr)->matrix[(index)])
#define RFUNCTION_SET(rf_ptr, index, value) ((rf_ptr)->matrix[(index)] = (value))
#define RFUNCTION_ADD(rf_ptr, index, value) ((rf_ptr)->matrix[(index)] += (value))

#define RFUNCTION_PRESUM_GET(rf_ptr, index)         ((rf_ptr)->presum_matrix[(index)])
#define RFUNCTION_PRESUM_SET(rf_ptr, index, value)  ((rf_ptr)->presum_matrix[(index)] = (value))
#define RFUNCTION_PRESUM_ADD(rf_ptr, index, value)  ((rf_ptr)->presum_matrix[(index)] += (value))

#define RFUNCTION_PEELNODE(rf_ptr)  ((rf_ptr)->peel_node)
#define RFUNCTION_TYPE(rf_ptr)      ((rf_ptr)->peel_type)

#define MAP_THETA(map_ptr, n)          ((map_ptr)->thetas[(n)])
#define MAP_INVERSETHETA(map_ptr, n)   ((map_ptr)->inversethetas[(n)])
#define MAP_HALF_THETA(map_ptr, n)          ((map_ptr)->halfthetas[(n)])
#define MAP_HALF_INVERSETHETA(map_ptr, n)   ((map_ptr)->halfinversethetas[(n)])
#define MAP_LENGTH(map_ptr)            ((map_ptr)->map_length)
#define MAP_PROB(map_ptr, n, val)      ((map_ptr)->markerprobs[((n) * 4) + val])
#define MAP_MINOR(map_ptr, n)          ((map_ptr)->allelefreqs[(n)])
#define MAP_MAJOR(map_ptr, n)          (1.0 - (map_ptr)->allelefreqs[(n)])

#define PERSON_ID(person_ptr)               ((person_ptr)->id)
#define PERSON_ISFOUNDER(person_ptr)        ((person_ptr)->isfounder)
#define PERSON_ISTYPED(person_ptr)          ((person_ptr)->istyped)
#define PERSON_DISEASEPROB(person_ptr, n)   ((person_ptr)->prob[(n)])
#define PERSON_GENOTYPE(person_ptr, n)      ((person_ptr)->genotypes[(n)])
#define PERSON_MOTHER(person_ptr)           ((person_ptr)->mother)
#define PERSON_FATHER(person_ptr)           ((person_ptr)->father)
#define PERSON_ISPARENT(person_ptr, n)      (((person_ptr)->mother == (n)) || ((person_ptr)->father == (n)))

#define DESCENTGRAPH_OFFSET(dg_ptr, personid, locusid, parentid) \
    (((dg_ptr)->subgraph_length * (locusid)) + ((personid) * 2) + (parentid))

#define DESCENTGRAPH_GET(dg_ptr, n) \
    ((dg_ptr)->graph[(n)])

#define DESCENTGRAPH_SET(dg_ptr, n, value) \
    ((dg_ptr)->graph[(n)] = (value))


#ifdef __cplusplus
extern "C" { 
#endif    
    void run_gpu_print_kernel(struct gpu_state* state);
    void run_gpu_print_pedigree_kernel(struct person* p, int length);
    
    void run_gpu_lsampler_kernel(int numblocks, int numthreads, struct gpu_state* state, int window_length, int offset);
    void run_gpu_lsampler_onepeel_kernel(int numblocks, int numthreads, struct gpu_state* state, int offset, int function_index);
    void run_gpu_lsampler_sample_kernel(int numblocks, int numthreads, struct gpu_state* state, int offset);
    
    void run_gpu_msampler_kernel(int numblocks, int numthreads, struct gpu_state* state, int meiosis);
    void run_gpu_msampler_reset_kernel(int numblocks, int numthreads, struct gpu_state* state, int meiosis, size_t shared);
    void run_gpu_msampler_likelihood_kernel(int numblocks, int numthreads, struct gpu_state* state, int meiosis, int last_meiosis, size_t shared);
    void run_gpu_msampler_sampling_kernel(struct gpu_state* state, int meiosis);
    void run_gpu_msampler_window_sampling_kernel(int numblocks, int numthreads, struct gpu_state* state, int meiosis, int offset);
    
    void run_gpu_lodscore_kernel(int numblocks, int numthreads, struct gpu_state* state);
    void run_gpu_lodscore_onepeel_kernel(int numblocks, int numthreads, struct gpu_state* state, int function_offset);
    void run_gpu_lodscoreinit_kernel(int numblocks, double* lodscores);
    void run_gpu_lodscorenormalise_kernel(int numblocks, struct gpu_state* state, int count, double trait_likelihood);
    void run_gpu_lodscoreprint_kernel(struct gpu_state* state);
    
    void run_gpu_curand_init_kernel(int numblocks, curandState* states, long int* seeds);
    void run_gpu_tinymt_init_kernel(int numblocks, tinymt32_status_t* states, uint32_t* params, uint32_t* seeds);
    
    void setup_lsampler_kernel();
    void setup_lodscore_kernel();
    void setup_msampler_kernel();
#ifdef __cplusplus
}
#endif


#endif

