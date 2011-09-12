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
#include "person.h"
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
    size_t mem_per_sampler;
    size_t mem_pedigree;
    size_t mem_map;
    vector<PeelOperation>& ops = psg.get_peel_order();
    
    mem_map = sizeof(struct geneticmap) + \
                (2 * sizeof(float) * (map->num_markers() - 1)) + \
                (4 * sizeof(float) *  map->num_markers());
    
    mem_pedigree = ped->num_members() * sizeof(struct person);
    for(unsigned i = 0; i < ped->num_members(); ++i) {
        Person* tmp = ped->get_by_index(i);
        
        if(tmp->istyped()) {
            mem_pedigree += (sizeof(int) * tmp->num_markers());
        }
    }
    
    mem_per_sampler = 0;
    for(unsigned i = 0; i < ops.size(); ++i) {
        double s = ops[i].get_cutset_size();
        
        mem_per_sampler += (pow(4.0, s)     * sizeof(float));   // matrix
        mem_per_sampler += (pow(4.0, s + 1) * sizeof(float));   // presum_matrix
        mem_per_sampler += (int(s + 1)      * sizeof(int));     // cutset
        mem_per_sampler += sizeof(struct rfunction);            // containing struct
    }
    
    return (num_samplers() * mem_per_sampler) + mem_pedigree + mem_map + sizeof(struct gpu_state);
}

unsigned GPUWrapper::num_samplers() {
    //return (map->num_markers() + 1) / 2; // + 1 in case there is an odd number of markers in the map
    return map->num_markers(); // wasteful yes, but let's keep everything easy to address ;-P
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

void GPUWrapper::kill_everything() {
    // rfunctions
    for(int i = 0; i < data->functions_length; ++i) {
        struct rfunction* rf = &(data->functions[i]);
        
        free(rf->cutset);
        free(rf->matrix);
        free(rf->presum_matrix);
    }
    
    free(data->functions);
    
    // pedigree
    for(int i = 0; i < data->pedigree_length; ++i) {
        struct person* p = &(data->pedigree[i]);
        
        free(p->genotypes);
    }
    
    free(data->pedigree);
    
    // map info
    free(data->map->thetas);
    free(data->map->inversethetas);
    free(data->map->markerprobs);
    free(data->map);
    
    // descent graph copy
    free(data->dg->graph);
    free(data->dg);
    
    // state
    free(data);
}

void GPUWrapper::init(PeelSequenceGenerator& psg) {
    
    size_t mem_needed = calculate_memory_requirements(psg);
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
    
    data = (struct gpu_state*) malloc(sizeof(struct gpu_state));
    if(!data) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    init_rfunctions(psg);
    init_pedigree();
    init_map();
    init_descentgraph();
}

void GPUWrapper::init_descentgraph() {
    
    data->dg = (struct descentgraph*) malloc(sizeof(struct descentgraph));
    if(!(data->dg)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    data->dg->subgraph_length = 2 * ped->num_members();
    data->dg->graph_length = 2 * ped->num_members() * map->num_markers();
    
    data->dg->graph = (int*) malloc(sizeof(int) * (2 * ped->num_members() * map->num_markers()));
    if(!(data->dg->graph)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
}

void GPUWrapper::init_map() {
    
    data->map = (struct geneticmap*) malloc(sizeof(struct geneticmap));
    if(!(data->map)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }

    data->map->map_length = map->num_markers();
    
    data->map->thetas = (float*) malloc(sizeof(float) * (map->num_markers() - 1));
    if(!(data->map->thetas)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    data->map->inversethetas = (float*) malloc(sizeof(float) * (map->num_markers() - 1));
    if(!(data->map->inversethetas)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    data->map->markerprobs = (float*) malloc(sizeof(float) * 4 * map->num_markers());
    if(!(data->map->markerprobs)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    for(unsigned i = 0; i < (map->num_markers() - 1); ++i) {
        data->map->thetas[i] = map->get_theta(i);
        data->map->inversethetas[i] = map->get_inversetheta(i);
    }
    
    for(unsigned i = 0; i < map->num_markers(); ++i) {
        Snp& marker = map->get_marker(i);
        for(unsigned j = 0; j < 4; ++j) {
            enum phased_trait pt = static_cast<enum phased_trait>(j);
            data->map->markerprobs[(i * 4) + j] = marker.get_prob(pt);
        }
    }
}

void GPUWrapper::init_pedigree() {

    data->pedigree_length = ped->num_members();
    data->pedigree = (struct person*) malloc(ped->num_members() * sizeof(struct person));
    if(!(data->pedigree)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    for(unsigned i = 0; i < ped->num_members(); ++i) {
        Person* tmp = ped->get_by_index(i);
        struct person* p = &(data->pedigree[i]);
        
        p->id = tmp->get_internalid();
        p->mother = tmp->get_maternalid();
        p->father = tmp->get_paternalid();
        p->isfounder = tmp->isfounder() ? 1 : 0;
        p->istyped = tmp->istyped() ? 1 : 0;
        p->genotypes = NULL;
        
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

void GPUWrapper::init_rfunctions(PeelSequenceGenerator& psg) {
    
    vector<PeelOperation>& ops = psg.get_peel_order();
    
    // need to allocate a shed-load of r-functions
    data->functions_length = num_samplers() * ops.size();
    data->functions_per_locus = ops.size();
    data->functions = (struct rfunction*) malloc(num_samplers() * ops.size() * sizeof(struct rfunction));
    if(!data->functions) {
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
        
        int prev1_index = -1;
        int prev2_index = -1;
        
        find_previous_functions(ops, j, prev1_index, prev2_index);
        
        if(prev1_index >= int(j)) {
            fprintf(stderr, "error: rfunction %d : prev1 %d\n", j, prev1_index);
            abort();
        }
        
        if(prev2_index >= int(j)) {
            fprintf(stderr, "error: rfunction %d : prev1 %d\n", j, prev2_index);
            abort();
        }
        
        for(unsigned i = 0; i < num_samp; ++i) {
            struct rfunction* rf = &data->functions[(i * num_func_per_samp) + j];
            
            rf->peel_type = peel_type;
            rf->peel_node = peel_node;
            rf->prev1 = (prev1_index == -1) ? NULL : &data->functions[(i * num_func_per_samp) + prev1_index];
            rf->prev2 = (prev2_index == -1) ? NULL : &data->functions[(i * num_func_per_samp) + prev2_index];
            
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

// -----

void GPUWrapper::find_previous_functions(vector<PeelOperation>& ops, int current_index, int& prev1_index, int& prev2_index) {
    if(ops[current_index].get_type() == CHILD_PEEL) {
        find_child_functions(ops, current_index, prev1_index, prev2_index);
    }
    else {
        find_generic_functions(ops, current_index, prev1_index, prev2_index);
    }
}

void GPUWrapper::find_generic_functions(vector<PeelOperation>& ops, int current_index, int& prev1_index, int& prev2_index) {
    vector<unsigned> tmp;
    tmp.push_back(ops[current_index].get_peelnode());
    
    prev1_index = find_function_containing(ops, current_index, tmp);
    prev2_index = find_function_containing(ops, current_index, tmp);
}

void GPUWrapper::find_child_functions(vector<PeelOperation>& ops, int current_index, int& prev1_index, int& prev2_index) {
    Person* p = ped->get_by_index(ops[current_index].get_peelnode());
    vector<unsigned> tmp;
    tmp.push_back(p->get_maternalid());
    tmp.push_back(p->get_paternalid());
    
    prev1_index = find_function_containing(ops, current_index, tmp);
    
    // don't even bother looking if the child is a leaf
    if(p->isleaf()) {
        prev2_index = -1;
        return;
    }
    
    tmp.clear();
    tmp.push_back(ops[current_index].get_peelnode());
    
    prev2_index = find_function_containing(ops, current_index, tmp);
}

int GPUWrapper::find_function_containing(vector<PeelOperation>& ops, int current_index, vector<unsigned>& nodes) {
    for(int i = 0; i < int(ops.size()); ++i) {
        if(ops[i].is_used()) {
            continue;
        }
        
        if(i >= current_index) {
            break;
        }
        
        if(ops[i].contains_cutnodes(nodes)) {
            ops[i].set_used();
            return i;
        }
    }
    
    return -1;
}

void GPUWrapper::copy_to_gpu(DescentGraph& dg) {
    for(unsigned i = 0; i < unsigned(data->dg->graph_length); ++i) {
        data->dg->graph[i] = dg.get_bit(i);
    }
}

void GPUWrapper::copy_from_gpu(DescentGraph& dg) {
    for(unsigned i = 0; i < unsigned(data->dg->graph_length); ++i) {
        dg.set_bit(i, data->dg->graph[i]);
    }
}

void GPUWrapper::step(DescentGraph& dg) {
    copy_to_gpu(dg);
    
    for(unsigned i = 0; i < map->num_markers(); ++i) {
        sampler_step(data, i);
    }
    
    copy_from_gpu(dg);
}

