using namespace std;

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <cmath>

#include "gpu_rfunction.h"
#include "gpu_wrapper.h"

#include "peeling.h"
#include "pedigree.h"
#include "genetic_map.h"
#include "peel_sequence_generator.h"
#include "descent_graph.h"


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
size_t GPUWrapper::calculate_memory_requirements(PeelSequenceGenerator& psg) {
    size_t mem_per_sampler = 0;
    vector<PeelOperation>& ops = psg.get_peel_order();
    
    for(unsigned i = 0; i < ops.size(); ++i) {
        double s = ops[i].get_cutset_size();
        
        mem_per_sampler += (pow(4.0, s)     * sizeof(float));   // matrix
        mem_per_sampler += (pow(4.0, s + 1) * sizeof(float));   // presum_matrix
        mem_per_sampler += (int(s + 1)      * sizeof(int));     // cutset
        mem_per_sampler += sizeof(struct rfunction);            // containing struct
    }
    
    return (num_samplers() * mem_per_sampler);
}

unsigned GPUWrapper::num_samplers() {
    return (map->num_markers() + 1) / 2; // + 1 in case there is an odd number of markers in the map
}

int GPUWrapper::convert_type(enum peeloperation type) {
    
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

void GPUWrapper::init(PeelSequenceGenerator& psg) {
    
    vector<PeelOperation>& ops = psg.get_peel_order();
    
    size_t mem_needed = calculate_memory_requirements(psg);
    /*
    if(mem_needed > SOME_DISCOVERED_LIMIT) {
        fprintf(stderr, 
                "error: your GPU does not have enough memory to hold everything! (needed: %d MB, found: %d MB)\n",
                mem_needed, SOME_DISCOVERED_LIMIT);
        abort();
    }
    */
    
    fprintf(stderr, "%d MB required for GPU rfunctions (%d bytes)\n", int(mem_needed / 1e6), int(mem_needed));
    return;
    
    // need to allocate a shed-load of r-functions
    data = (struct rfunction*) malloc(num_samplers() * ops.size() * sizeof(struct rfunction));
    if(!data) {
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
    
    for(unsigned j = 0; j < ops.size(); ++j) {
        double s = ops[j].get_cutset_size();
        int peel_type = convert_type(ops[j].get_type());
        int peel_node = ops[j].get_peelnode();
        int matrix_length = static_cast<int>(pow(4.0, s));
        int presum_length = static_cast<int>(pow(4.0, s + 1));
        int cutset_length = ops[j].get_cutset_size() + 1;
        
        for(unsigned i = 0; i < num_samp; ++i) {
            struct rfunction* rf = &data[(i * num_samp) + j];
            
            rf->peel_type = peel_type;
            rf->peel_node = peel_node;
            
            rf->matrix_length = matrix_length;
            rf->matrix = (float*) malloc(matrix_length * sizeof(float));
            if(! rf->matrix) {
                fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
                abort();
            }
            
            rf->presum_length = presum_length;
            rf->presum_matrix = (float*) malloc(presum_length * sizeof(float));
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
        }
    }
    
}

