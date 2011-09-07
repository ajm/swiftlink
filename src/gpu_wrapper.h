#ifndef LKG_GPUWRAPPER_H_
#define LKG_GPUWRAPPER_H_

#include <cstdlib>

#include "peeling.h"
#include "gpu_rfunction.h"

class Pedigree;
class GeneticMap;
class PeelSequenceGenerator;
class DescentGraph;

class GPUWrapper {
    
    Pedigree* ped;
    GeneticMap* map;
    
    struct gpu_state* data;
    
    size_t calculate_memory_requirements(PeelSequenceGenerator& psg);
    unsigned num_samplers();
    int convert_type(enum peeloperation type);
    void init(PeelSequenceGenerator& psg);
    void init_map();
    void init_pedigree();
    void init_rfunctions(PeelSequenceGenerator& psg);
    
    void kill_everything();
    
    void find_previous_functions(vector<PeelOperation>& ops, int current_index, int& prev1_index, int& prev2_index);
    void find_generic_functions(vector<PeelOperation>& ops, int current_index, int& prev1_index, int& prev2_index);
    void find_child_functions(vector<PeelOperation>& ops, int current_index, int& prev1_index, int& prev2_index);
    int find_function_containing(vector<PeelOperation>& ops, int current_index, vector<unsigned>& nodes);
    
 public :
    GPUWrapper(Pedigree* ped, GeneticMap* map, PeelSequenceGenerator& psg) :
        ped(ped),
        map(map),
        data(NULL) {
        
        init(psg);    
    }
        
    GPUWrapper(const GPUWrapper& rhs) :
        ped(rhs.ped),
        map(rhs.map),
        data(rhs.data) {}
    
    ~GPUWrapper() {
        kill_everything();
    }
    
    GPUWrapper& operator=(const GPUWrapper& rhs) {
        
        if(&rhs != this) {
            ped = rhs.ped;
            map = rhs.map;
            data = rhs.data; // XXX this could be a problem, but I don't keep PeelSequenceGenerator obj
        }
        
        return *this;
    }
    
    void step(DescentGraph& dg, unsigned parameter) {}
};

#endif

