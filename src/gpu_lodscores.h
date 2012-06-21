#ifndef LKG_GPUWRAPPER_H_
#define LKG_GPUWRAPPER_H_

#include <cstdlib>

#include "cuda_quiet.h"

#include "peel_sequence_generator.h"
#include "peeling.h"
#include "types.h"


class Pedigree;
class GeneticMap;
class PeelSequenceGenerator;
class DescentGraph;
class LODscores;

class GPULodscores {
    
    Pedigree* ped;
    GeneticMap* map;
    PeelSequenceGenerator* psg;
    struct mcmc_options options;
    
    struct gpu_state* loc_state;
    struct gpu_state* dev_state;
    int* dev_graph;
    double* dev_lodscores;
    
    unsigned int num_lodscore_threads;
    unsigned int count;
    double trait_likelihood;
    
    int lodscore_num_blocks();
    int convert_type(enum peeloperation type);
    
    void init();
    void init_map();
    void init_pedigree();
    void init_descentgraph();
    void init_founderallelegraph();
    
    void gpu_init(vector<PeelOperation>& ops);
    struct geneticmap* gpu_init_map();
    struct person* gpu_init_pedigree();
    struct rfunction* gpu_init_rfunctions(vector<PeelOperation>& ops);
    struct descentgraph* gpu_init_descentgraph();
    double* gpu_init_lodscores();
    
    
    void copy_to_gpu(DescentGraph& dg);
    void copy_from_gpu(DescentGraph& dg);
    
    void select_best_gpu();
    
    void kill_everything();
    
    int optimal_lodscore_threads();
    
    
 public :
    GPULodscores(Pedigree* ped, GeneticMap* map, PeelSequenceGenerator* psg, struct mcmc_options options, double trait_prob) :
        ped(ped),
        map(map),
        psg(psg),
        options(options),
        loc_state(NULL),
        dev_state(NULL),
        dev_graph(NULL),
        dev_lodscores(NULL),
        num_lodscore_threads(128),
        count(0),
        trait_likelihood(trait_prob) {
        
        vector<PeelOperation>& ops = psg->get_peel_order();
        
        init();
        gpu_init(ops);
    }
        
    GPULodscores(const GPULodscores& rhs) :
        ped(rhs.ped),
        map(rhs.map),
        psg(rhs.psg),
        options(rhs.options),
        loc_state(rhs.loc_state),
        dev_state(rhs.dev_state),
        dev_graph(rhs.dev_graph),
        dev_lodscores(rhs.dev_lodscores),
        num_lodscore_threads(rhs.num_lodscore_threads),
        count(rhs.count),
        trait_likelihood(rhs.trait_likelihood) {}
    
    ~GPULodscores() {
        kill_everything();
        
        cudaDeviceReset();
    }
    
    GPULodscores& operator=(const GPULodscores& rhs) {
        if(&rhs != this) {
            ped = rhs.ped;
            map = rhs.map;
            psg = rhs.psg;
            options = rhs.options;
            loc_state = rhs.loc_state;
            dev_state = rhs.dev_state;
            dev_graph = rhs.dev_graph;
            dev_lodscores = rhs.dev_lodscores;
            num_lodscore_threads = rhs.num_lodscore_threads;
            count = rhs.count;
            trait_likelihood = rhs.trait_likelihood;
        }
        return *this;
    }
    
    void calculate(DescentGraph& dg);
    void block_until_finished();
    void get_results(LODscores* lod);
};

#endif

