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
    
    struct rfunction* data;
    
    size_t calculate_memory_requirements(PeelSequenceGenerator& psg);
    unsigned num_samplers();
    int convert_type(enum peeloperation type);
    void init(PeelSequenceGenerator& psg);
    
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
        // XXX
    }
    
    GPUWrapper& operator=(const GPUWrapper& rhs) {
        
        if(&rhs != this) {
            ped = rhs.ped;
            map = rhs.map;
            data = rhs.data;
        }
        
        return *this;
    }
    
    void step(DescentGraph& dg, unsigned parameter) {}
};

#endif

