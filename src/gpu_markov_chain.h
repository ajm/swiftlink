#ifndef LKG_GPUWRAPPER_H_
#define LKG_GPUWRAPPER_H_

#include <cstdlib>

#include "cuda_quiet.h"

#include "peeling.h"
#include "peel_sequence_generator.h"
#include "types.h"
//#include "tinymt/tinymt32_host.h"


#define fp_type double

class Pedigree;
class GeneticMap;
class PeelSequenceGenerator;
class DescentGraph;
class LODscores;

class GPUMarkovChain {
    
    Pedigree* ped;
    GeneticMap* map;
    PeelSequenceGenerator* psg;
    struct mcmc_options options;
    
    struct gpu_state* loc_state;
    struct gpu_state* dev_state;
    int* dev_graph;
    fp_type* dev_lodscores;
    
    size_t calculate_memory_requirements(vector<PeelOperation>& ops);
    unsigned num_samplers();
    int num_threads_per_block();
    int num_blocks();
    int lsampler_num_blocks();
    int msampler_num_blocks();
    int lodscore_num_blocks();
    int windowed_msampler_blocks(int window_length);
    int convert_type(enum peeloperation type);
    void find_founderallelegraph_ordering(struct gpu_state* state);
    
    void init(vector<PeelOperation>& ops);
    void init_map();
    void init_pedigree();
    void init_rfunctions(vector<PeelOperation>& ops);
    void init_descentgraph();
    void init_founderallelegraph();
    
    void gpu_init(vector<PeelOperation>& ops);
    struct geneticmap* gpu_init_map();
    struct person* gpu_init_pedigree();
    struct rfunction* gpu_init_rfunctions(vector<PeelOperation>& ops);
    struct descentgraph* gpu_init_descentgraph();
    curandState* gpu_init_random_curand();
    //tinymt32_status_t* gpu_init_random_tinymt();
    double* gpu_init_lodscores();
    struct founderallelegraph* gpu_init_founderallelegraph();
    
    void print_person(struct person* p);
    
    void copy_to_gpu(DescentGraph& dg);
    void copy_from_gpu(DescentGraph& dg);
    
    void select_best_gpu();
    
    void kill_everything();
    
    int optimal_lodscore_threads();
    int optimal_lsampler_threads();
    
    
 public :
    GPUMarkovChain(Pedigree* ped, GeneticMap* map, PeelSequenceGenerator* psg, struct mcmc_options options) :
        ped(ped),
        map(map),
        psg(psg),
        options(options),
        loc_state(NULL),
        dev_state(NULL),
        dev_graph(NULL),
        dev_lodscores(NULL) {
        
        vector<PeelOperation>& ops = psg->get_peel_order();
        
        init(ops);
        gpu_init(ops);
        
        //print_everything(loc_state);
        //run_gpu_print_kernel(dev_state);
    }
        
    GPUMarkovChain(const GPUMarkovChain& rhs) :
        ped(rhs.ped),
        map(rhs.map),
        psg(rhs.psg),
        options(rhs.options),
        loc_state(rhs.loc_state),
        dev_state(rhs.dev_state),
        dev_graph(rhs.dev_graph),
        dev_lodscores(rhs.dev_lodscores) {}
    
    ~GPUMarkovChain() {
        kill_everything();
        
        cudaDeviceReset();
    }
    
    GPUMarkovChain& operator=(const GPUMarkovChain& rhs) {
        if(&rhs != this) {
            ped = rhs.ped;
            map = rhs.map;
            psg = rhs.psg;
            options = rhs.options;
            loc_state = rhs.loc_state; // XXX this could be a problem, but I don't keep PeelSequenceGenerator obj
            dev_state = rhs.dev_state;
            dev_graph = rhs.dev_graph;
            dev_lodscores = rhs.dev_lodscores;
        }
        return *this;
    }
    
    LODscores* run(DescentGraph& dg);
};

#endif

