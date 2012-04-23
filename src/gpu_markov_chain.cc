using namespace std;

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <cmath>
#include <ctime>
#include <cfloat>

#include <iostream>
#include <fstream>

#include "cuda_quiet.h"
#include "cuda_common.h"
//#include "tinymt/tinymt32_host.h"

#include "gpu_markov_chain.h"

#include "peeling.h"
#include "pedigree.h"
#include "person.h"
#include "genetic_map.h"
#include "descent_graph.h"
#include "progress.h"
#include "random.h"

#include "peeler.h"
#include <vector>

#define CUDA_CALLANDTEST(x) \
do {\
    if((x) != cudaSuccess) {\
        fprintf(stderr, "error: %s (%s:%d %s())\n", cudaGetErrorString(cudaGetLastError()), __FILE__, __LINE__, __func__);\
        cudaDeviceReset();\
        abort();\
    }\
} while(0)


/*
 * TODO
 * 
 * I need a routine that will calculate the total amount of memory
 * required for this, then (in the proper GPU version) compare that
 * to what is actually available
 *
 * At the moment if there is not enough memory just abort() and bring
 * up a firefox instance pointing to amazon or nvidia or something...
 *
 * however it would be *really* cool to just work out the limit then 
 * run the maximum number concurrently and have the code just automatically
 * do the rest once they are finished (would be amazingly cool as MCMCMC
 * would become easier to implement)
 */
size_t GPUMarkovChain::calculate_memory_requirements(vector<PeelOperation>& ops) {
    size_t mem_per_sampler;
    size_t mem_pedigree;
    size_t mem_map;
    size_t mem_rand;
    size_t mem_fag;
    
    mem_map = sizeof(struct geneticmap) + \
                (4 * sizeof(fp_type) * (map->num_markers() - 1)) + \
                (4 * sizeof(fp_type) *  map->num_markers());
    
    mem_pedigree = ped->num_members() * sizeof(struct person);
    for(unsigned i = 0; i < ped->num_members(); ++i) {
        Person* tmp = ped->get_by_index(i);
        
        if(tmp->istyped()) {
            mem_pedigree += (sizeof(int) * tmp->num_markers());
        }
    }
    
    mem_per_sampler = 0;
    //mem_shared = 0;
    for(unsigned i = 0; i < ops.size(); ++i) {
        double s = ops[i].get_cutset_size();
        
        mem_per_sampler += (pow(4.0, s)     * sizeof(fp_type));   // matrix
        mem_per_sampler += (pow(4.0, s + 1) * sizeof(fp_type));   // presum_matrix
        
        mem_per_sampler += (int(s + 1)      * sizeof(int));     // cutset
        mem_per_sampler += sizeof(struct rfunction);            // containing struct
    }
    
    mem_rand = sizeof(curandState) * num_blocks();
    
    int num_fa = ped->num_founders() * 2;
    //mem_fag = (sizeof(struct founderallelegraph) + (num_fa * sizeof(int)) + (num_fa * num_fa * sizeof(struct adjacent_node))) * map->num_markers();
    
    return (num_samplers() * mem_per_sampler) + mem_pedigree + mem_map + mem_rand + mem_fag + sizeof(struct gpu_state) + (sizeof(int) * ped->num_members());
}

unsigned GPUMarkovChain::num_samplers() {
    return map->num_markers();
}

int GPUMarkovChain::num_threads_per_block() {
    //return NUM_THREADS;
    return 96;
}

int GPUMarkovChain::lsampler_num_blocks() {
    return (map->num_markers() / 2) + ((map->num_markers() % 2) == 0 ? 0 : 1);
}

int GPUMarkovChain::num_blocks() {
    return lsampler_num_blocks();
}

int GPUMarkovChain::msampler_num_blocks() {
    return (map->num_markers() / 8) + ((map->num_markers() % 8) == 0 ? 0 : 1);
    // 8 because I will use 256 threads
    // 256 (threads) / 32 (threads per warp) = 8
}

int GPUMarkovChain::lodscore_num_blocks() {
    return map->num_markers() - 1;
}

int GPUMarkovChain::convert_type(enum peeloperation type) {
    
    switch(type) {
        case CHILD_PEEL:
            return GPU_CHILD_PEEL;
        case PARENT_PEEL:
            return GPU_PARENT_PEEL;
        case PARTNER_PEEL:
        case LAST_PEEL:
            return GPU_PARTNER_PEEL;
        default :
            break;
    }
    
    abort();
}

void GPUMarkovChain::kill_everything() {
    // rfunctions
    for(int i = 0; i < loc_state->functions_length; ++i) {
        struct rfunction* rf = &(loc_state->functions[i]);
        
        free(rf->cutset);
        free(rf->matrix);
        free(rf->presum_matrix);
    }
    
    free(loc_state->functions);
    
    // pedigree
    for(int i = 0; i < loc_state->pedigree_length; ++i) {
        struct person* p = &(loc_state->pedigree[i]);
        
        free(p->genotypes);
    }
    
    free(loc_state->pedigree);
    
    // map info
    free(loc_state->map->thetas);
    free(loc_state->map->inversethetas);
    free(loc_state->map->markerprobs);
    free(loc_state->map);
    
    // descent graph copy
    free(loc_state->dg->graph);
    free(loc_state->dg);
    
    // state
    free(loc_state);
}

void GPUMarkovChain::init(vector<PeelOperation>& ops) {
    
    size_t mem_needed = calculate_memory_requirements(ops);
    /*
    if(mem_needed > SOME_DISCOVERED_LIMIT) {
        fprintf(stderr, 
                "error: your GPU does not have enough memory to hold everything! (needed: %d MB, found: %d MB)\n",
                mem_needed, SOME_DISCOVERED_LIMIT);
        abort();
    }
    */
    
    fprintf(stderr, "%.2f MB required for GPU rfunctions (%d bytes)\n", mem_needed / 1e6, int(mem_needed));
    //return;
    
    loc_state = (struct gpu_state*) malloc(sizeof(struct gpu_state));
    if(!loc_state) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    init_rfunctions(ops);
    init_pedigree();
    init_map();
    init_descentgraph();
    init_founderallelegraph();
}

void GPUMarkovChain::gpu_init(vector<PeelOperation>& ops) {
    select_best_gpu();
    
    struct gpu_state tmp;
    
    tmp.map = gpu_init_map(); fprintf(stderr, "GPU: map loaded\n");
    tmp.dg = gpu_init_descentgraph(); fprintf(stderr, "GPU: descent graph loaded\n");
    
    tmp.pedigree = gpu_init_pedigree(); fprintf(stderr, "GPU: pedigree loaded\n");
    tmp.pedigree_length = loc_state->pedigree_length;
    
    tmp.functions = gpu_init_rfunctions(ops); fprintf(stderr, "GPU: rfunctions loaded\n");
    tmp.functions_length = loc_state->functions_length;
    tmp.functions_per_locus = loc_state->functions_per_locus;
    
    tmp.randstates = gpu_init_random_curand(); fprintf(stderr, "GPU: random loaded\n");
    //tmp.randmt = gpu_init_random_tinymt();
    
    tmp.lodscores = gpu_init_lodscores(); fprintf(stderr, "GPU: lod scores loaded\n");
    
    //tmp.graphs = gpu_init_founderallelegraph();
    tmp.founderallele_count = loc_state->founderallele_count;
    
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp.raw_matrix, sizeof(double) * map->num_markers() * 2));
    fprintf(stderr, "GPU: msampler matrix loaded\n");
    
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp.fa_sequence, sizeof(int) * ped->num_members()));
    CUDA_CALLANDTEST(cudaMemcpy(tmp.fa_sequence, loc_state->fa_sequence, sizeof(int) * ped->num_members(), cudaMemcpyHostToDevice));
    fprintf(stderr, "GPU: msampler order loaded\n");
    
    CUDA_CALLANDTEST(cudaMalloc((void**)&dev_state, sizeof(struct gpu_state)));
    CUDA_CALLANDTEST(cudaMemcpy(dev_state, &tmp, sizeof(struct gpu_state), cudaMemcpyHostToDevice));
    fprintf(stderr, "GPU: memory loaded\n");
}

curandState* GPUMarkovChain::gpu_init_random_curand() {
    long int* host_seeds;
    long int* dev_seeds;
    curandState* dev_rng_state;
    int prng_num = num_blocks();
    fstream randfile;
    
    host_seeds = (long int*) malloc(sizeof(*host_seeds) * prng_num);
    if(! host_seeds) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    randfile.open ("/dev/urandom", ios::in);
    if(not randfile.is_open()) {
        fprintf(stderr, "error: could not open /dev/urandom\n");
        abort();
    }
    
    //printf("sizeof(*host_seeds) = %d\n", sizeof(host_seeds));
    
    for(int i = 0; i < prng_num; ++i) {
        //seeds[i] = get_random();
        long int data;
        
        randfile.read((char*)&data, sizeof(data));
        
        host_seeds[i] = data;
    }
    
    randfile.close();
    
    CUDA_CALLANDTEST(cudaMalloc((void**) &dev_seeds,   sizeof(*dev_seeds) * prng_num));
    CUDA_CALLANDTEST(cudaMemcpy(dev_seeds, host_seeds, sizeof(*dev_seeds) * prng_num, cudaMemcpyHostToDevice));
    
    CUDA_CALLANDTEST(cudaMalloc((void**) &dev_rng_state, sizeof(*dev_rng_state) * prng_num));
    
    
    run_gpu_curand_init_kernel(num_blocks(), dev_rng_state, dev_seeds);
    
    
    free(host_seeds);
    cudaFree(dev_seeds);
    
    return dev_rng_state;
}
/*
tinymt32_status_t* GPUMarkovChain::gpu_init_random_tinymt() {
    uint32_t* params;
    uint32_t* dparams;
    uint32_t* seeds;
    uint32_t* dseeds;
    tinymt32_status_t* dev_rng_state;
    int prng_num = num_blocks();
    fstream randfile;
    
    params = (uint32_t*) malloc(sizeof(*params) * 3 * prng_num);
    if(! params) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    // this file was generated by TinyMT32DC
    #define TINYMT_PARAM_FILENAME "tinymt32dc.0.65536.txt"
    
    if(tinymt32_set_params(TINYMT_PARAM_FILENAME, params, prng_num) != 0) {
        fprintf(stderr, "error: could not read in parameters from %s\n", TINYMT_PARAM_FILENAME);
        abort();
    }
    
    seeds = (uint32_t*) malloc(sizeof(*seeds) * prng_num);
    if(! seeds) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    randfile.open ("/dev/urandom", ios::in);
    if(not randfile.is_open()) {
        fprintf(stderr, "error: could not open /dev/random\n");
        abort();
    }
    
    for(int i = 0; i < prng_num; ++i) {
        //seeds[i] = get_random();
        uint32_t data;
        
        randfile.read((char*)&data, 4);
        
        seeds[i] = data;
    }
    
    randfile.close();
    
    CUDA_CALLANDTEST(cudaMalloc((void**) &dparams, sizeof(*dparams) * 3 * prng_num));
    CUDA_CALLANDTEST(cudaMemcpy(dparams, params, sizeof(*dparams) * 3 * prng_num, cudaMemcpyHostToDevice));
    
    CUDA_CALLANDTEST(cudaMalloc((void**) &dseeds, sizeof(*dseeds) * prng_num));
    CUDA_CALLANDTEST(cudaMemcpy(dseeds, seeds, sizeof(*dseeds) * prng_num, cudaMemcpyHostToDevice));
    
    
    CUDA_CALLANDTEST(cudaMalloc((void**) &dev_rng_state, sizeof(*dev_rng_state) * prng_num));
    
    
    run_gpu_tinymt_init_kernel(num_blocks(), dev_rng_state, dparams, dseeds);
    
    cudaThreadSynchronize();
    
    
    free(params);
    free(seeds);
    cudaFree(dparams);
    cudaFree(dseeds);
    
    #undef TINYMT_PARAM_FILENAME
    
    return dev_rng_state;
}
*/
void GPUMarkovChain::init_founderallelegraph() {

    int num_founderalleles = ped->num_founders() * 2;
    /*
    loc_state->graphs = (struct founderallelegraph*) malloc(sizeof(struct founderallelegraph) * map->num_markers());
    if(!(loc_state->graphs)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    */
    
    /*
    for(unsigned i = 0; i < map->num_markers(); ++i) {
        loc_state->graphs[i].num_neighbours = (int*) malloc(sizeof(int) * num_founderalleles);
        loc_state->graphs[i].graph = (struct adjacent_node*) malloc(sizeof(struct adjacent_node*) * num_founderalleles * num_founderalleles);        
    }
    */
    
    loc_state->founderallele_count = num_founderalleles;
    
    loc_state->fa_sequence = (int*) malloc(sizeof(int) * ped->num_members());
    if(!(loc_state->fa_sequence)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    find_founderallelegraph_ordering(loc_state);
}

/*
struct founderallelegraph* GPUMarkovChain::gpu_init_founderallelegraph() {
    
    //int num_founderalleles = ped->num_founders() * 2;
    
    struct founderallelegraph* tmp = (struct founderallelegraph*) malloc(sizeof(struct founderallelegraph) * map->num_markers());
    if(!tmp) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    
    
    //for(unsigned i = 0; i < map->num_markers(); ++i) {
    //    CUDA_CALLANDTEST(cudaMalloc((void**)&(tmp[i].num_neighbours), sizeof(int) * num_founderalleles));
    //    CUDA_CALLANDTEST(cudaMalloc((void**)&(tmp[i].graph), sizeof(struct adjacent_node) * num_founderalleles * num_founderalleles));
    //}
    
    
    struct founderallelegraph* dev_fag;
    
    CUDA_CALLANDTEST(cudaMalloc((void**) &dev_fag, sizeof(struct founderallelegraph) * map->num_markers()));
    CUDA_CALLANDTEST(cudaMemcpy(dev_fag, tmp,    sizeof(struct founderallelegraph) * map->num_markers(), cudaMemcpyHostToDevice));
    
    free(tmp);
    
    return dev_fag;
}
*/
void GPUMarkovChain::init_descentgraph() {
    
    loc_state->dg = (struct descentgraph*) malloc(sizeof(struct descentgraph));
    if(!(loc_state->dg)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    loc_state->dg->subgraph_length = 2 * ped->num_members();
    loc_state->dg->graph_length = 2 * ped->num_members() * map->num_markers();
    loc_state->dg->transmission_prob = log(0.5) * (2 * (ped->num_members() - ped->num_founders()));
    
    loc_state->dg->graph = (int*) malloc(sizeof(int) * (2 * ped->num_members() * map->num_markers()));
    if(!(loc_state->dg->graph)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
}

struct descentgraph* GPUMarkovChain::gpu_init_descentgraph() {

    struct descentgraph tmp;
    struct descentgraph* dev_dg;
    
    tmp.graph_length = loc_state->dg->graph_length;
    tmp.subgraph_length = loc_state->dg->subgraph_length;
    tmp.transmission_prob = loc_state->dg->transmission_prob;
    
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp.graph, sizeof(int) * tmp.graph_length));
    CUDA_CALLANDTEST(cudaMemcpy(tmp.graph, loc_state->dg->graph, sizeof(int) * tmp.graph_length, cudaMemcpyHostToDevice));
    
    CUDA_CALLANDTEST(cudaMalloc((void**)&dev_dg, sizeof(struct descentgraph)));
    CUDA_CALLANDTEST(cudaMemcpy(dev_dg, &tmp, sizeof(struct descentgraph), cudaMemcpyHostToDevice));
    
    dev_graph = tmp.graph;
    
    return dev_dg;
}

void GPUMarkovChain::init_map() {
    
    loc_state->map = (struct geneticmap*) malloc(sizeof(struct geneticmap));
    if(!(loc_state->map)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }

    loc_state->map->map_length = map->num_markers();
    
    loc_state->map->thetas = (fp_type*) malloc(sizeof(fp_type) * (map->num_markers() - 1));
    if(!(loc_state->map->thetas)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    loc_state->map->inversethetas = (fp_type*) malloc(sizeof(fp_type) * (map->num_markers() - 1));
    if(!(loc_state->map->inversethetas)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    loc_state->map->halfthetas = (fp_type*) malloc(sizeof(fp_type) * (map->num_markers() - 1));
    if(!(loc_state->map->halfthetas)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    loc_state->map->halfinversethetas = (fp_type*) malloc(sizeof(fp_type) * (map->num_markers() - 1));
    if(!(loc_state->map->halfinversethetas)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    loc_state->map->markerprobs = (fp_type*) malloc(sizeof(fp_type) * 4 * map->num_markers());
    if(!(loc_state->map->markerprobs)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }

    loc_state->map->allelefreqs = (fp_type*) malloc(sizeof(fp_type) * map->num_markers());
    if(!(loc_state->map->allelefreqs)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    for(unsigned i = 0; i < (map->num_markers() - 1); ++i) {
        loc_state->map->thetas[i] = map->get_theta(i);
        loc_state->map->inversethetas[i] = map->get_inversetheta(i);
        loc_state->map->halfthetas[i] = map->get_theta_halfway(i);
        loc_state->map->halfinversethetas[i] = map->get_inversetheta_halfway(i);
    }
    
    for(unsigned i = 0; i < map->num_markers(); ++i) {
        Snp& marker = map->get_marker(i);
        for(unsigned j = 0; j < 4; ++j) {
            enum phased_trait pt = static_cast<enum phased_trait>(j);
            loc_state->map->markerprobs[(i * 4) + j] = marker.get_prob(pt);
        }
        loc_state->map->allelefreqs[i] = marker.minor();
    }
}

struct geneticmap* GPUMarkovChain::gpu_init_map() {
    
    struct geneticmap tmp;
    struct geneticmap* dev_map;
    
    tmp.map_length = loc_state->map->map_length;
    
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp.thetas,            sizeof(fp_type) * (map->num_markers() - 1)));
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp.inversethetas,     sizeof(fp_type) * (map->num_markers() - 1)));
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp.halfthetas,        sizeof(fp_type) * (map->num_markers() - 1)));
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp.halfinversethetas, sizeof(fp_type) * (map->num_markers() - 1)));
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp.markerprobs,       sizeof(fp_type) * (map->num_markers() * 4)));
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp.allelefreqs,       sizeof(fp_type) *  map->num_markers()));
    
    CUDA_CALLANDTEST(cudaMemcpy(tmp.thetas,             loc_state->map->thetas,             sizeof(fp_type) * (map->num_markers() - 1), cudaMemcpyHostToDevice));
    CUDA_CALLANDTEST(cudaMemcpy(tmp.inversethetas,      loc_state->map->inversethetas,      sizeof(fp_type) * (map->num_markers() - 1), cudaMemcpyHostToDevice));
    CUDA_CALLANDTEST(cudaMemcpy(tmp.halfthetas,         loc_state->map->halfthetas,         sizeof(fp_type) * (map->num_markers() - 1), cudaMemcpyHostToDevice));
    CUDA_CALLANDTEST(cudaMemcpy(tmp.halfinversethetas,  loc_state->map->halfinversethetas,  sizeof(fp_type) * (map->num_markers() - 1), cudaMemcpyHostToDevice));
    CUDA_CALLANDTEST(cudaMemcpy(tmp.markerprobs,        loc_state->map->markerprobs,        sizeof(fp_type) * (map->num_markers() * 4), cudaMemcpyHostToDevice));
    CUDA_CALLANDTEST(cudaMemcpy(tmp.allelefreqs,        loc_state->map->allelefreqs,        sizeof(fp_type) *  map->num_markers(),      cudaMemcpyHostToDevice));
    
    CUDA_CALLANDTEST(cudaMalloc((void**)&dev_map, sizeof(struct geneticmap)));
    CUDA_CALLANDTEST(cudaMemcpy(dev_map, &tmp, sizeof(struct geneticmap), cudaMemcpyHostToDevice));
    
    return dev_map;
}

void GPUMarkovChain::init_pedigree() {

    loc_state->pedigree_length = ped->num_members();
    loc_state->pedigree = (struct person*) malloc(ped->num_members() * sizeof(struct person));
    if(!(loc_state->pedigree)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    for(unsigned i = 0; i < ped->num_members(); ++i) {
        Person* tmp = ped->get_by_index(i);
        struct person* p = &(loc_state->pedigree[i]);
        
        p->id = tmp->get_internalid();
        p->mother = tmp->get_maternalid();
        p->father = tmp->get_paternalid();
        p->isfounder = tmp->isfounder() ? 1 : 0;
        p->istyped = tmp->istyped() ? 1 : 0;
        p->genotypes = NULL;
        p->genotypes_length = tmp->istyped() ? map->num_markers() : 0 ;
        
        // probs
        for(int j = 0; j < 4; ++j) {
            p->prob[j] = tmp->get_disease_prob(static_cast<enum phased_trait>(j));
        }
        
        // genotypes
        if(tmp->istyped()) {
            p->genotypes = (int*) malloc(ped->num_markers() * sizeof(int));
            if(!(p->genotypes)) {
                fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
                abort();
            }
            
            for(unsigned j = 0; j < ped->num_markers(); ++j) {
                p->genotypes[j] = static_cast<int>(tmp->get_marker(j));
            }
        }
    }
}

struct person* GPUMarkovChain::gpu_init_pedigree() {
    
    struct person tmp;
    struct person* dev_ped;
    
    CUDA_CALLANDTEST(cudaMalloc((void**) &dev_ped, sizeof(struct person) * loc_state->pedigree_length));
    
    /*
    printf("host: sizeof(struct person) = %d\n", sizeof(struct person));
    
    run_gpu_print_pedigree_kernel(dev_ped, loc_state->pedigree_length); // XXX
    cudaThreadSynchronize();
    printf("====================\n");
    */
    
    for(int i = 0; i < loc_state->pedigree_length; ++i) {
        struct person* p = &(loc_state->pedigree[i]);
        struct person* dev_p = &(dev_ped[i]);
        
        memcpy(&tmp, p, sizeof(struct person));
        
        //print_person(&tmp);
        
        CUDA_CALLANDTEST(cudaMalloc((void**) &(tmp.genotypes), sizeof(int) * p->genotypes_length));
        CUDA_CALLANDTEST(cudaMemcpy(tmp.genotypes, p->genotypes, sizeof(int) * p->genotypes_length, cudaMemcpyHostToDevice));
        
        CUDA_CALLANDTEST(cudaMemcpy(dev_p, &tmp, sizeof(struct person), cudaMemcpyHostToDevice));
    }
    
    /*
    printf("====================\n");
    run_gpu_print_pedigree_kernel(dev_ped, loc_state->pedigree_length); // XXX
    cudaThreadSynchronize();
    */
    
    return dev_ped;
}

void GPUMarkovChain::init_rfunctions(vector<PeelOperation>& ops) {
    
    // need to allocate a shed-load of r-functions
    loc_state->functions_length = num_samplers() * ops.size();
    loc_state->functions_per_locus = ops.size();
    loc_state->functions = (struct rfunction*) malloc(num_samplers() * ops.size() * sizeof(struct rfunction));
    if(!loc_state->functions) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    // allocate the memory for each r-function
    //  - matrix + matrix_length
    //  - presum_matrix + presum_length
    //  - cutset + cutset_length
    // fill in the details for each r-function
    //  - cutset
    //  - peel_node
    //  - peel_type
    //  - prev1 + prev2 (more complex, but only worked out once, then just repeated in all with offsets recalculated)
    unsigned num_samp = num_samplers();
    unsigned num_func_per_samp = ops.size();
    
    for(unsigned j = 0; j < num_func_per_samp; ++j) {
        double s = ops[j].get_cutset_size();
        int peel_type = convert_type(ops[j].get_type());
        int peel_node = ops[j].get_peelnode();
        int matrix_length = static_cast<int>(pow(4.0, s));
        int presum_length = static_cast<int>(pow(4.0, s + 1));
        int cutset_length = ops[j].get_cutset_size() + 1;
        
        int prev1_index = ops[j].get_previous_op1();
        int prev2_index = ops[j].get_previous_op2();
        
        for(unsigned i = 0; i < num_samp; ++i) {
            vector<int>* vindices = ops[j].get_valid_indices(i);
        
            struct rfunction* rf = &loc_state->functions[(i * num_func_per_samp) + j];
            
            rf->id = j;
            rf->peel_type = peel_type;
            rf->peel_node = peel_node;
            rf->prev1 = (prev1_index == -1) ? NULL : &loc_state->functions[(i * num_func_per_samp) + prev1_index];
            rf->prev2 = (prev2_index == -1) ? NULL : &loc_state->functions[(i * num_func_per_samp) + prev2_index];
            
            rf->matrix_length = matrix_length;
            rf->matrix = (fp_type*) malloc(matrix_length * sizeof(fp_type));
            if(! rf->matrix) {
                fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
                abort();
            }
            
            rf->presum_length = presum_length;
            rf->presum_matrix = (fp_type*) malloc(presum_length * sizeof(fp_type));
            if(! rf->presum_matrix) {
                fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
                abort();
            }
            
            rf->cutset_length = cutset_length;
            rf->cutset = (int*) malloc(cutset_length * sizeof(int));
            if(! rf->cutset) {
                fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
                abort();
            }
            
            for(unsigned k = 0; k < unsigned(cutset_length - 1); ++k) {
                rf->cutset[k] = ops[j].get_cutnode(k);
            }
            rf->cutset[cutset_length - 1] = peel_node;
            
            // valid indices from genotype elimination
            rf->presum_indices_length = vindices->size();
            rf->presum_indices = (int*) malloc(rf->presum_indices_length * sizeof(int));
            if(! rf->presum_indices) {
                fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
                abort();
            }
            
            copy(vindices->begin(), vindices->end(), rf->presum_indices);
        }
    }   
}

struct rfunction* GPUMarkovChain::gpu_init_rfunctions(vector<PeelOperation>& ops) {
    
    struct rfunction tmp;
    struct rfunction* dev_rfunc;
    
    CUDA_CALLANDTEST(cudaMalloc((void**) &dev_rfunc, sizeof(struct rfunction) * ops.size() * map->num_markers()));
    
    
    for(int i = 0; i < int(ops.size()); ++i) {
        // get all the stuff for this rfunction
        int prev1_index = ops[i].get_previous_op1();
        int prev2_index = ops[i].get_previous_op2();
        
        int children_length = ops[i].get_children_size();
        int* children = NULL;
        
        if(children_length > 0) {
            children = new int[children_length];
            copy(ops[i].get_children().begin(), ops[i].get_children().end(), children);
        }
        
        int cutset_length = ops[i].get_cutset_size() + 1;
        int* cutset = new int[cutset_length];
        for(int j = 0; j < (cutset_length - 1); ++j) {
            cutset[j] = ops[i].get_cutnode(j);
        }
        cutset[cutset_length-1] = ops[i].get_peelnode();
        
        tmp.id = i;
        tmp.peel_type = convert_type(ops[i].get_type());
        tmp.peel_node = ops[i].get_peelnode();
        
        tmp.cutset_length = cutset_length;
        tmp.presum_length = int(pow(4.0, ops[i].get_cutset_size() + 1));
        tmp.matrix_length = int(pow(4.0, ops[i].get_cutset_size()));
        tmp.children_length = children_length;
        
        
        fprintf(stderr, "  rfunctions: loading %d...\n", i);
        
        for(int j = 0; j < int(map->num_markers()); ++j) {
            // 3 - 5 mallocs, 1 - 2 memcpys
            CUDA_CALLANDTEST(cudaMalloc((void**) &tmp.presum_matrix, sizeof(fp_type) * tmp.presum_length));
            CUDA_CALLANDTEST(cudaMemset(tmp.presum_matrix, 0, sizeof(fp_type) * tmp.presum_length));
            
            CUDA_CALLANDTEST(cudaMalloc((void**) &tmp.matrix,        sizeof(fp_type) * tmp.matrix_length));
            CUDA_CALLANDTEST(cudaMalloc((void**) &tmp.cutset,        sizeof(int)     * tmp.cutset_length));
            CUDA_CALLANDTEST(cudaMemcpy(tmp.cutset, cutset, sizeof(int) * tmp.cutset_length, cudaMemcpyHostToDevice));
            
            if(children_length > 0) {
                CUDA_CALLANDTEST(cudaMalloc((void**) &tmp.children, sizeof(int) * tmp.children_length));
                CUDA_CALLANDTEST(cudaMemcpy(tmp.children, children, sizeof(int) * tmp.children_length, cudaMemcpyHostToDevice));
            
                CUDA_CALLANDTEST(cudaMalloc((void**) &tmp.transmission, sizeof(fp_type) * tmp.children_length * 64));
            }
            else {
                tmp.children = NULL;
                tmp.transmission = NULL;
            }
            
            // 2 pointers
            tmp.prev1 = (prev1_index != -1) ? &dev_rfunc[(j * ops.size()) + prev1_index] : NULL;
            tmp.prev2 = (prev2_index != -1) ? &dev_rfunc[(j * ops.size()) + prev2_index] : NULL;
            
            
            // valid indices from genotype elimination
            vector<int>* vindices = ops[i].get_valid_indices(j);
            
            tmp.presum_indices_length = vindices->size();
            //fprintf(stderr, "length = %d valid = %d\n", tmp.presum_length, tmp.valid_indices_length);
            int* presum_indices = new int[vindices->size()];
            copy(vindices->begin(), vindices->end(), presum_indices);
            
            CUDA_CALLANDTEST(cudaMalloc((void**) &tmp.presum_indices, sizeof(int) * tmp.presum_indices_length));
            CUDA_CALLANDTEST(cudaMemcpy(tmp.presum_indices, presum_indices, sizeof(int) * tmp.presum_indices_length, cudaMemcpyHostToDevice));
            
            delete[] presum_indices;
            // 
            
            // copy tmp to device via dev_rf pointer
            struct rfunction* dev_rf = &dev_rfunc[(j * ops.size()) + i];
            
            CUDA_CALLANDTEST(cudaMemcpy(dev_rf, &tmp, sizeof(struct rfunction), cudaMemcpyHostToDevice));
        }
        
        if(children_length > 0)
            delete[] children;
        
        delete[] cutset;
    }
    
    return dev_rfunc;
}

fp_type* GPUMarkovChain::gpu_init_lodscores() {
    fp_type* tmp;
    
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp, sizeof(fp_type) * (map->num_markers() - 1)));
    
    run_gpu_lodscoreinit_kernel(map->num_markers() - 1, tmp);
    
    cudaThreadSynchronize();
    
    dev_lodscores = tmp;
    
    return tmp;
}

void GPUMarkovChain::copy_to_gpu(DescentGraph& dg) {
    for(unsigned i = 0; i < unsigned(loc_state->dg->graph_length); ++i) {
        loc_state->dg->graph[i] = dg.get_bit(i);
    }
}

void GPUMarkovChain::copy_from_gpu(DescentGraph& dg) {
    for(unsigned i = 0; i < unsigned(loc_state->dg->graph_length); ++i) {
        dg.set_bit(i, loc_state->dg->graph[i]);
    }
}

int GPUMarkovChain::windowed_msampler_blocks(int window_length) {
    return (map->num_markers() / window_length) + ((map->num_markers() % window_length) == 0 ? 0 : 1);
}

int GPUMarkovChain::optimal_lodscore_threads() {
        
    int threadcount[] = {32, 64, 96, 128, 160, 192, 224, 256};
    
    float mintime = FLT_MAX;
    int minthread = 32;
    
    for(int i = 0; i < 8; ++i) {
    
        float total = 0.0;
    
        for(int j = 0; j < 10; ++j) {
            cudaEvent_t start;
            cudaEvent_t stop;
            
            cudaEventCreate(&start);
            cudaEventCreate(&stop);
            
            cudaEventRecord(start, 0);
            
            
            run_gpu_lodscore_kernel(lodscore_num_blocks(), threadcount[i], dev_state);
            cudaThreadSynchronize();
            
            cudaError_t error = cudaGetLastError();
            if(error != cudaSuccess) {
                printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
                abort();
            }
            
            
            cudaEventRecord(stop, 0);
            cudaEventSynchronize(stop);
            
            float elapsedTime;
            cudaEventElapsedTime(&elapsedTime, start, stop);
            
            total += elapsedTime;
            
            cudaEventDestroy(start);
            cudaEventDestroy(stop);
        }
        
        total /= 10.0;
        
        fprintf(stderr, "lodscores: %d threads - %.3fms\n", threadcount[i], total);
        
        if(total < mintime) {
            mintime = total;
            minthread = threadcount[i];
        }
    }
    
    fprintf(stderr, "optimal LOD score threads = %d\n", minthread);
    
    return minthread;
}

int GPUMarkovChain::optimal_lsampler_threads() {
    int threadcount[] = {32, 64, 96, 128, 160, 192, 224, 256};
    
    float mintime = FLT_MAX;
    int minthread = 32;
    
    vector<int> l_ordering;
    
    int markers_per_window = 2;
    for(int i = 0; i < markers_per_window; ++i)
        l_ordering.push_back(i);
    
    
    for(int i = 0; i < 8; ++i) {
    
        float total = 0.0;
    
        for(int j = 0; j < 10; ++j) {
            cudaEvent_t start;
            cudaEvent_t stop;
            
            cudaEventCreate(&start);
            cudaEventCreate(&stop);
            
            cudaEventRecord(start, 0);
            
            
            for(int k = 0; k < int(l_ordering.size()); ++k) {
                int num_blocks = (map->num_markers() / markers_per_window) + \
                                ((map->num_markers() % markers_per_window == 0) ? 0 : 1);
                
                run_gpu_lsampler_kernel(num_blocks, threadcount[i], dev_state, markers_per_window, l_ordering[k]);
                cudaThreadSynchronize();
                
                cudaError_t error = cudaGetLastError();
                if(error != cudaSuccess) {
                    printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
                    abort();
                }
            }
            
            
            cudaEventRecord(stop, 0);
            cudaEventSynchronize(stop);
            
            float elapsedTime;
            cudaEventElapsedTime(&elapsedTime, start, stop);
            
            total += elapsedTime;
            
            cudaEventDestroy(start);
            cudaEventDestroy(stop);
        }
        
        total /= 10.0;
        
        fprintf(stderr, "lsampler: %d threads - %.3fms\n", threadcount[i], total);
        
        if(total < mintime) {
            mintime = total;
            minthread = threadcount[i];
        }
    }
    
    fprintf(stderr, "optimal locus sampler threads = %d\n", minthread);
    
    return minthread;
}


double* GPUMarkovChain::run(DescentGraph& dg, double trait_likelihood) {
    int count = 0;
    int even_count;
    int odd_count;
    cudaError_t error;
    int num_meioses = 2 * (ped->num_members() - ped->num_founders());
    //size_t shared_mem = 8 * ((2 * (2 * ped->num_founders()) * (2 * ped->num_founders())) + (2 * ped->num_founders())) ;
    size_t shared_mem = (14 * (2 * ped->num_founders())) + (2 * ped->num_members());
    shared_mem += (shared_mem % 4);
    shared_mem *= 8;
    
    //int msampler_blocks = (map->num_markers() / 10) + ((map->num_markers() % 10) == 0 ? 0 : 1);
    
    //int window_index = 0;
    //int window_length[] = {21, 23, 25};
    int window_length = 20;
    //int window_length = 32;
    
    
    //Peeler peel(ped, map, psg, 0);
    //double trait_likelihood = peel.get_trait_prob();
    
    vector<PeelOperation>& ops = psg->get_peel_order();
    
    /*
    run_gpu_print_kernel(dev_state);
    cudaThreadSynchronize();
    error = cudaGetLastError();
    if(error != cudaSuccess) {
        printf("CUDA kernel error: %s\n", cudaGetErrorString(error));
        abort();
    }
    exit(0);
    */
    
    odd_count  = map->num_markers() / 2;
    even_count = odd_count + ((map->num_markers() % 2) != 0 ? 1 : 0);
    
    setup_lsampler_kernel();
    setup_lodscore_kernel();
    setup_msampler_kernel();
    
    fprintf(stderr, "GPU: setup done\n");
    
    
    //run_gpu_print_kernel(dev_state);
    
    
    
    
    CUDA_CALLANDTEST(cudaMemcpy(dev_graph, dg.get_internal_ptr(), dg.get_internal_size(), cudaMemcpyHostToDevice));
    
    int lodscore_threads = optimal_lodscore_threads();
    int lsampler_threads = optimal_lsampler_threads();
    
    //abort();
    
    
    vector<int> l_ordering;
    vector<int> m_ordering;
    
    int markers_per_window = 2;
    for(int i = 0; i < markers_per_window; ++i)
        l_ordering.push_back(i);
    
    //for(int i = 0; i < num_meioses; ++i)
    //    m_ordering.push_back(i);
        
    for(unsigned int i = 0; i < num_meioses; ++i) {
        unsigned person_id = ped->num_founders() + (i / 2);
        enum parentage p = static_cast<enum parentage>(i % 2);
        
        Person* tmp = ped->get_by_index(person_id);
        
        if(not tmp->safe_to_ignore_meiosis(p)) {
            m_ordering.push_back(i);
        }
    }
    
        
    printf("requesting %.2f KB (%d bytes) of shared memory per block\n", shared_mem / 1024.0, (int)shared_mem);
    
    CUDA_CALLANDTEST(cudaMemcpy(dev_graph, dg.get_internal_ptr(), dg.get_internal_size(), cudaMemcpyHostToDevice));
    
    Progress p("CUDA MCMC: ", options.iterations + options.burnin);
    
    
    for(int i = 0; i < (options.iterations + options.burnin); ++i) {
        
        if(get_random() < options.lsampler_prob) {
        //if((i % 2) == 0) {
            
            /*
            run_gpu_lsampler_kernel(even_count, num_threads_per_block(), dev_state, 2, 0);            
            cudaThreadSynchronize();
            error = cudaGetLastError();
            if(error != cudaSuccess) {
                printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
                abort();
            }
            
            
            run_gpu_lsampler_kernel(odd_count, num_threads_per_block(), dev_state, 2, 1);
            cudaThreadSynchronize();
            error = cudaGetLastError();
            if(error != cudaSuccess) {
                printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
                abort();
            }
            */
            
            
            random_shuffle(l_ordering.begin(), l_ordering.end());
            
//#define GPU_LSAMPLER_MICROBENCHMARKING 1
#ifdef GPU_LSAMPLER_MICROBENCHMARKING
            cudaEvent_t start;
            cudaEvent_t stop;
            
            cudaEventCreate(&start);
            cudaEventCreate(&stop);
            
            cudaEventRecord(start, 0);
#endif
            
            for(int j = 0; j < int(l_ordering.size()); ++j) {
                int num_blocks = (map->num_markers() / markers_per_window) + \
                                ((map->num_markers() % markers_per_window == 0) ? 0 : 1);
                
                run_gpu_lsampler_kernel(num_blocks, lsampler_threads, dev_state, markers_per_window, l_ordering[j]);
                
                cudaThreadSynchronize();
                error = cudaGetLastError();
                if(error != cudaSuccess) {
                    printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
                    abort();
                }
            }
            
#ifdef GPU_LSAMPLER_MICROBENCHMARKING
            cudaEventRecord(stop, 0);
            cudaEventSynchronize(stop);
            
            float elapsedTime; // in milliseconds
            cudaEventElapsedTime(&elapsedTime, start, stop);
            
            if(i < options.burnin)
                fprintf(stderr, "LSAMPLER %.3f\n", elapsedTime);
            
            cudaEventDestroy(start);
            cudaEventDestroy(stop);
#endif
            
            
            /*
            for(unsigned int j = 0; j < ops.size(); ++j) {
                int num_threads = static_cast<int>(pow(4.0, ops[j].get_cutset_size() + 1));
                
                run_gpu_lsampler_onepeel_kernel(even_count, num_threads, dev_state, 0, j);            
                cudaThreadSynchronize();
                error = cudaGetLastError();
                if(error != cudaSuccess) {
                    printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
                    abort();
                }
            }
            
            run_gpu_lsampler_sample_kernel(even_count, 32, dev_state, 0);
            cudaThreadSynchronize();
            error = cudaGetLastError();
            if(error != cudaSuccess) {
                printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
                abort();
            }
            
            
            for(unsigned int j = 0; j < ops.size(); ++j) {
                int num_threads = static_cast<int>(pow(4.0, ops[j].get_cutset_size() + 1));
                
                run_gpu_lsampler_onepeel_kernel(odd_count, num_threads, dev_state, 1, j);            
                cudaThreadSynchronize();
                error = cudaGetLastError();
                if(error != cudaSuccess) {
                    printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
                    abort();
                }
            }
            
            run_gpu_lsampler_sample_kernel(odd_count, 32, dev_state, 1);
            cudaThreadSynchronize();
            error = cudaGetLastError();
            if(error != cudaSuccess) {
                printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
                abort();
            }
            */
            
            
            
            /*
            // test legality of descent graph
            CUDA_CALLANDTEST(cudaMemcpy(dg.get_internal_ptr(), dev_graph, dg.get_internal_size(), cudaMemcpyDeviceToHost));
            
            if(dg.get_likelihood() == LOG_ILLEGAL) {
                fprintf(stderr, "error: descent graph illegal after l-sampler (%d)\n", i);
                abort();
            }
            */
        }
        else {
//#define GPU_MSAMPLER_MICROBENCHMARKING 1
#ifdef GPU_MSAMPLER_MICROBENCHMARKING
            cudaEvent_t start;
            cudaEvent_t stop;
            
            cudaEventCreate(&start);
            cudaEventCreate(&stop);
            
            cudaEventRecord(start, 0);
#endif

            random_shuffle(m_ordering.begin(), m_ordering.end());
            
            run_gpu_msampler_reset_kernel(msampler_num_blocks(), 256, dev_state, m_ordering[0], shared_mem);
            cudaThreadSynchronize();
            error = cudaGetLastError();
            if(error != cudaSuccess) {
                printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
                abort();
            }
            
            for(int j = 0; j < int(m_ordering.size()); ++j) {
                run_gpu_msampler_likelihood_kernel(msampler_num_blocks(), 256, dev_state, \
                                m_ordering[j], (j == 0) ? m_ordering[0] : m_ordering[j-1], shared_mem);
                cudaThreadSynchronize();
                error = cudaGetLastError();
                if(error != cudaSuccess) {
                    printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
                    abort();
                }
                
                //run_gpu_msampler_sampling_kernel(dev_state, m_ordering[j]);
                //run_gpu_msampler_window_sampling_kernel(windowed_msampler_blocks(window_length[window_index]), 32, dev_state, m_ordering[j], window_length[window_index]);
                run_gpu_msampler_window_sampling_kernel(windowed_msampler_blocks(window_length), 32, dev_state, m_ordering[j], window_length);
                cudaThreadSynchronize();
                error = cudaGetLastError();
                if(error != cudaSuccess) {
                    printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
                    abort();
                }
                
                //window_index = (window_index + 1) % 3;
                
                
                
                /*
                // test legality of descent graph
                CUDA_CALLANDTEST(cudaMemcpy(dg.get_internal_ptr(), dev_graph, dg.get_internal_size(), cudaMemcpyDeviceToHost));
                
                if(dg.get_likelihood() == LOG_ILLEGAL) {
                    fprintf(stderr, "error: descent graph illegal after m-sampler (meiosis %d) (%d)\n", j, i);
                    abort();
                }
                */
            }
#ifdef GPU_MSAMPLER_MICROBENCHMARKING
            cudaEventRecord(stop, 0);
            cudaEventSynchronize(stop);
            
            float elapsedTime; // in milliseconds
            cudaEventElapsedTime(&elapsedTime, start, stop);
            
            if(i < options.burnin)
                fprintf(stderr, "MSAMPLER %.3f\n", elapsedTime);
            
            cudaEventDestroy(start);
            cudaEventDestroy(stop);
#endif
        }
        
        
        p.increment();
        
        if(i < options.burnin)
            continue;
        
        if((i % options.scoring_period) == 0) {
            /*
            for(int i = 0; i < int(map->num_markers() - 1); ++i) {
                peelers[i]->process(&dg);
            }
            */
            
#ifdef GPU_LODSCORE_MICROBENCHMARKING
            cudaEvent_t start;
            cudaEvent_t stop;
            
            cudaEventCreate(&start);
            cudaEventCreate(&stop);
            
            cudaEventRecord(start, 0);
#endif
            
            /*
            for(unsigned int j = 0; j < ops.size(); ++j) {
                int num_threads = static_cast<int>(pow(4.0, ops[j].get_cutset_size() + 1));
                
                run_gpu_lodscore_onepeel_kernel(lodscore_num_blocks(), num_threads, dev_state, j);            
                cudaThreadSynchronize();
                error = cudaGetLastError();
                if(error != cudaSuccess) {
                    printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
                    abort();
                }
            }
            */
            
            run_gpu_lodscore_kernel(lodscore_num_blocks(), lodscore_threads, dev_state);
            cudaThreadSynchronize();
            
            cudaError_t error = cudaGetLastError();
            if(error != cudaSuccess) {
                printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
                abort();
            }
            
#ifdef GPU_LODSCORE_MICROBENCHMARKING
            cudaEventRecord(stop, 0);
            cudaEventSynchronize(stop);
            
            float elapsedTime; // in milliseconds
            cudaEventElapsedTime(&elapsedTime, start, stop);
            
            fprintf(stderr, "LODSCORE %.3f\n", elapsedTime);
            
            cudaEventDestroy(start);
            cudaEventDestroy(stop);
#endif
            
            count++;
            
        }
    }
    
    p.finish();
    
    
    run_gpu_lodscorenormalise_kernel(map->num_markers() - 1, dev_state, count, trait_likelihood);
    cudaThreadSynchronize();
    
    
    //run_gpu_lodscoreprint_kernel(dev_state);
    //cudaThreadSynchronize();
    //printf("count = %d, P(T) = %.3f\n", count, trait_likelihood / log(10));
    
    
    /*
    //double* lod_scores = new double[map->num_markers() - 1];
    for(unsigned int i = 0; i < (map->num_markers() - 1); ++i) {
        //lod_scores[i] = peelers[i]->get();
        printf("rs%d\t%f\n", i+1, peelers[i]->get());
        delete peelers[i];
    }
    //return lod_scores;
    */
    
    double* data = new double[map->num_markers() - 1];
    
    CUDA_CALLANDTEST(cudaMemcpy(data, dev_lodscores, (map->num_markers() - 1) * sizeof(double), cudaMemcpyDeviceToHost));
    
    return data;
}

void GPUMarkovChain::select_best_gpu() {
    int num_devices;
    int i;

    cudaGetDeviceCount(&num_devices);

    if(num_devices > 1) {
        int max_multiprocessors = 0;
        int max_device = 0;
        
        for(i = 0; i < num_devices; ++i) {
            cudaDeviceProp properties;
            cudaGetDeviceProperties(&properties, i);
            if(max_multiprocessors < properties.multiProcessorCount) {
                max_multiprocessors = properties.multiProcessorCount;
                max_device = i;
            }
        }
        
        cudaSetDevice(max_device);
    }
}

void GPUMarkovChain::find_founderallelegraph_ordering(struct gpu_state* state) {
    vector<unsigned> seq;
    vector<bool> visited(ped->num_members(), false);
    int total = ped->num_members();
    
    // we can start by putting in all founders as there are clearly
    // no dependencies
    for(unsigned i = 0; i < ped->num_founders(); ++i) {
        seq.push_back(i);
        visited[i] = true;
        total--;
    }
    
    while(total > 0) {
        for(unsigned i = ped->num_founders(); i < ped->num_members(); ++i) {
            if(visited[i])
                continue;
        
            Person* p = ped->get_by_index(i);
            
            if(visited[p->get_maternalid()] and visited[p->get_paternalid()]) {
                seq.push_back(i);
                visited[i] = true;
                total--;
            }
        }
    }
    
    if(seq.size() != ped->num_members()) {
        fprintf(stderr, "Founder allele sequence generation failed\n");
        abort();
    }
    
    for(unsigned i = 0; i < ped->num_members(); ++i) {
        state->fa_sequence[i] = static_cast<int>(seq[i]);
    }
}

void GPUMarkovChain::print_person(struct person* p) {
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
    
    //printf("\tgenotypes: ");
    //print_ints(p->genotypes, p->genotypes_length);
    printf("\n");
}

