#ifndef LKG_GPUWRAPPER_H_
#define LKG_GPUWRAPPER_H_

#include <cstdlib>

#include "cuda.h"
#include "cuda_runtime.h"

#include "peeling.h"
#include "peel_sequence_generator.h"
#include "gpu_rfunction.h"

class Pedigree;
class GeneticMap;
class PeelSequenceGenerator;
class DescentGraph;

class GPUWrapper {
    
    Pedigree* ped;
    GeneticMap* map;
    
    struct gpu_state* loc_state;
    struct gpu_state* dev_state;
    
    size_t calculate_memory_requirements(vector<PeelOperation>& ops);
    unsigned num_samplers();
    int convert_type(enum peeloperation type);
    
    void init(vector<PeelOperation>& ops);
    void init_map();
    void init_pedigree();
    void init_rfunctions(vector<PeelOperation>& ops);
    void init_descentgraph();
    
    void gpu_init(vector<PeelOperation>& ops);
    struct geneticmap* gpu_init_map();
    struct person* gpu_init_pedigree();
    struct rfunction* gpu_init_rfunctions(vector<PeelOperation>& ops);
    struct descentgraph* gpu_init_descentgraph();
    
    void copy_to_gpu(DescentGraph& dg);
    void copy_from_gpu(DescentGraph& dg);
    
    void select_best_gpu();
    
    void kill_everything();
    
    
 public :
    GPUWrapper(Pedigree* ped, GeneticMap* map, PeelSequenceGenerator& psg) :
        ped(ped),
        map(map),
        loc_state(NULL),
        dev_state(NULL) {
        
        vector<PeelOperation>& ops = psg.get_peel_order();
        init(ops);
        gpu_init(ops);
    }
        
    GPUWrapper(const GPUWrapper& rhs) :
        ped(rhs.ped),
        map(rhs.map),
        loc_state(rhs.loc_state),
        dev_state(rhs.dev_state) {}
    
    ~GPUWrapper() {
        kill_everything();
        
        cudaDeviceReset();
    }
    
    GPUWrapper& operator=(const GPUWrapper& rhs) {
        
        if(&rhs != this) {
            ped = rhs.ped;
            map = rhs.map;
            loc_state = rhs.loc_state; // XXX this could be a problem, but I don't keep PeelSequenceGenerator obj
            dev_state = rhs.dev_state;
        }
        
        return *this;
    }
    
    void step(DescentGraph& dg);
};

#endif

