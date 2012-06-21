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

#include "gpu_lodscores.h"

#include "peeling.h"
#include "pedigree.h"
#include "person.h"
#include "genetic_map.h"
#include "descent_graph.h"

#include "peeler.h"
#include "lod_score.h"

#include <vector>

#define CUDA_CALLANDTEST(x) \
do {\
    if((x) != cudaSuccess) {\
        fprintf(stderr, "error: %s (%s:%d %s())\n", cudaGetErrorString(cudaGetLastError()), __FILE__, __LINE__, __func__);\
        cudaDeviceReset();\
        abort();\
    }\
} while(0)


void GPULodscores::select_best_gpu() {
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

int GPULodscores::convert_type(enum peeloperation type) {
    
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

void GPULodscores::kill_everything() {
    
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
    free(loc_state->map->partialthetas);
    free(loc_state->map);
    
    // descent graph copy
    free(loc_state->dg->graph);
    free(loc_state->dg);
    
    // state
    free(loc_state);
}

void GPULodscores::init() {
    
    loc_state = (struct gpu_state*) malloc(sizeof(struct gpu_state));
    if(!loc_state) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    init_pedigree();
    init_map();
    init_descentgraph();
}

void GPULodscores::gpu_init(vector<PeelOperation>& ops) {
    select_best_gpu();
    
    struct gpu_state tmp;
    
    tmp.map = gpu_init_map(); fprintf(stderr, "GPU: map loaded\n");
    tmp.dg = gpu_init_descentgraph(); fprintf(stderr, "GPU: descent graph loaded\n");
    
    tmp.pedigree = gpu_init_pedigree(); fprintf(stderr, "GPU: pedigree loaded\n");
    tmp.pedigree_length = loc_state->pedigree_length;
    
    tmp.functions = gpu_init_rfunctions(ops); fprintf(stderr, "GPU: rfunctions loaded\n");
    tmp.functions_length = map->num_markers() * ops.size();
    tmp.functions_per_locus = ops.size();
    
    tmp.lodscores = gpu_init_lodscores(); fprintf(stderr, "GPU: lod scores loaded\n");
    tmp.lodscores_per_marker = map->get_lodscore_count();
    
    CUDA_CALLANDTEST(cudaMalloc((void**)&dev_state, sizeof(struct gpu_state)));
    CUDA_CALLANDTEST(cudaMemcpy(dev_state, &tmp, sizeof(struct gpu_state), cudaMemcpyHostToDevice));
    fprintf(stderr, "GPU: memory loaded\n");
    
    
    
    setup_lodscore_kernel();
    num_lodscore_threads = optimal_lodscore_threads();
    
    run_gpu_lodscoreinit_kernel(map->num_markers() - 1, map->get_lodscore_count(), tmp.lodscores, map->get_lodscore_count());
    
    cudaThreadSynchronize();
    
    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess) {
        printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
        abort();
    }
    
    size_t avail, total;
    cudaMemGetInfo( &avail, &total );
    size_t used = total - avail;
    
    printf("Device memory used: %f MB\n", used / 1e6);
}

void GPULodscores::init_descentgraph() {
    
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

struct descentgraph* GPULodscores::gpu_init_descentgraph() {

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

void GPULodscores::init_map() {
    
    loc_state->map = (struct geneticmap*) malloc(sizeof(struct geneticmap));
    if(!(loc_state->map)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }

    loc_state->map->map_length = map->num_markers();
    loc_state->map->lodscores_per_marker = map->get_lodscore_count();
    
    loc_state->map->thetas = (double*) malloc(sizeof(double) * (map->num_markers() - 1));
    if(!(loc_state->map->thetas)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    loc_state->map->inversethetas = (double*) malloc(sizeof(double) * (map->num_markers() - 1));
    if(!(loc_state->map->inversethetas)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    loc_state->map->halfthetas = (double*) malloc(sizeof(double) * (map->num_markers() - 1));
    if(!(loc_state->map->halfthetas)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    loc_state->map->halfinversethetas = (double*) malloc(sizeof(double) * (map->num_markers() - 1));
    if(!(loc_state->map->halfinversethetas)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    loc_state->map->partialthetas = (double*) malloc(sizeof(double) * (map->num_markers() - 1));
    if(!(loc_state->map->partialthetas)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    loc_state->map->markerprobs = (double*) malloc(sizeof(double) * 4 * map->num_markers());
    if(!(loc_state->map->markerprobs)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }

    loc_state->map->allelefreqs = (double*) malloc(sizeof(double) * map->num_markers());
    if(!(loc_state->map->allelefreqs)) {
        fprintf(stderr, "error: %s (%s:%d)\n", strerror(errno), __FILE__, __LINE__);
        abort();
    }
    
    for(unsigned i = 0; i < (map->num_markers() - 1); ++i) {
        loc_state->map->thetas[i] = map->get_theta(i);
        loc_state->map->inversethetas[i] = map->get_inversetheta(i);
        loc_state->map->halfthetas[i] = map->get_theta_halfway(i);
        loc_state->map->halfinversethetas[i] = map->get_inversetheta_halfway(i);
        loc_state->map->partialthetas[i] = map->get_theta_partial_raw(i);
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

struct geneticmap* GPULodscores::gpu_init_map() {
    
    struct geneticmap tmp;
    struct geneticmap* dev_map;
    
    tmp.map_length = loc_state->map->map_length;
    tmp.lodscores_per_marker = loc_state->map->lodscores_per_marker;
    
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp.thetas,            sizeof(double) * (map->num_markers() - 1)));
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp.inversethetas,     sizeof(double) * (map->num_markers() - 1)));
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp.halfthetas,        sizeof(double) * (map->num_markers() - 1)));
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp.halfinversethetas, sizeof(double) * (map->num_markers() - 1)));
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp.partialthetas,     sizeof(double) * (map->num_markers() - 1)));
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp.markerprobs,       sizeof(double) * (map->num_markers() * 4)));
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp.allelefreqs,       sizeof(double) *  map->num_markers()));
    
    CUDA_CALLANDTEST(cudaMemcpy(tmp.thetas,             loc_state->map->thetas,             sizeof(double) * (map->num_markers() - 1), cudaMemcpyHostToDevice));
    CUDA_CALLANDTEST(cudaMemcpy(tmp.inversethetas,      loc_state->map->inversethetas,      sizeof(double) * (map->num_markers() - 1), cudaMemcpyHostToDevice));
    CUDA_CALLANDTEST(cudaMemcpy(tmp.halfthetas,         loc_state->map->halfthetas,         sizeof(double) * (map->num_markers() - 1), cudaMemcpyHostToDevice));
    CUDA_CALLANDTEST(cudaMemcpy(tmp.halfinversethetas,  loc_state->map->halfinversethetas,  sizeof(double) * (map->num_markers() - 1), cudaMemcpyHostToDevice));
    CUDA_CALLANDTEST(cudaMemcpy(tmp.partialthetas,      loc_state->map->partialthetas,      sizeof(double) * (map->num_markers() - 1), cudaMemcpyHostToDevice));
    CUDA_CALLANDTEST(cudaMemcpy(tmp.markerprobs,        loc_state->map->markerprobs,        sizeof(double) * (map->num_markers() * 4), cudaMemcpyHostToDevice));
    CUDA_CALLANDTEST(cudaMemcpy(tmp.allelefreqs,        loc_state->map->allelefreqs,        sizeof(double) *  map->num_markers(),      cudaMemcpyHostToDevice));
    
    CUDA_CALLANDTEST(cudaMalloc((void**)&dev_map, sizeof(struct geneticmap)));
    CUDA_CALLANDTEST(cudaMemcpy(dev_map, &tmp, sizeof(struct geneticmap), cudaMemcpyHostToDevice));
    
    return dev_map;
}

void GPULodscores::init_pedigree() {

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

struct person* GPULodscores::gpu_init_pedigree() {
    
    struct person tmp;
    struct person* dev_ped;
    
    CUDA_CALLANDTEST(cudaMalloc((void**) &dev_ped, sizeof(struct person) * loc_state->pedigree_length));
    
    for(int i = 0; i < loc_state->pedigree_length; ++i) {
        struct person* p = &(loc_state->pedigree[i]);
        struct person* dev_p = &(dev_ped[i]);
        
        memcpy(&tmp, p, sizeof(struct person));
        
        CUDA_CALLANDTEST(cudaMalloc((void**) &(tmp.genotypes), sizeof(int) * p->genotypes_length));
        CUDA_CALLANDTEST(cudaMemcpy(tmp.genotypes, p->genotypes, sizeof(int) * p->genotypes_length, cudaMemcpyHostToDevice));
        
        CUDA_CALLANDTEST(cudaMemcpy(dev_p, &tmp, sizeof(struct person), cudaMemcpyHostToDevice));
    }
    
    return dev_ped;
}

struct rfunction* GPULodscores::gpu_init_rfunctions(vector<PeelOperation>& ops) {
    
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
            CUDA_CALLANDTEST(cudaMalloc((void**) &tmp.presum_matrix, sizeof(double) * tmp.presum_length));
            CUDA_CALLANDTEST(cudaMemset(tmp.presum_matrix, 0, sizeof(double) * tmp.presum_length));
            
            CUDA_CALLANDTEST(cudaMalloc((void**) &tmp.matrix,        sizeof(double) * tmp.matrix_length));
            CUDA_CALLANDTEST(cudaMalloc((void**) &tmp.cutset,        sizeof(int)     * tmp.cutset_length));
            CUDA_CALLANDTEST(cudaMemcpy(tmp.cutset, cutset, sizeof(int) * tmp.cutset_length, cudaMemcpyHostToDevice));
            
            if(children_length > 0) {
                CUDA_CALLANDTEST(cudaMalloc((void**) &tmp.children, sizeof(int) * tmp.children_length));
                CUDA_CALLANDTEST(cudaMemcpy(tmp.children, children, sizeof(int) * tmp.children_length, cudaMemcpyHostToDevice));
            
                CUDA_CALLANDTEST(cudaMalloc((void**) &tmp.transmission, sizeof(double) * tmp.children_length * 64));
            }
            else {
                tmp.children = NULL;
                tmp.transmission = NULL;
            }
            
            // 2 pointers
            tmp.prev1 = (prev1_index != -1) ? &dev_rfunc[(j * ops.size()) + prev1_index] : NULL;
            tmp.prev2 = (prev2_index != -1) ? &dev_rfunc[(j * ops.size()) + prev2_index] : NULL;
            
            
            // valid indices from genotype elimination
            vector<int>* vindices = ops[i].get_presum_indices(j);
            
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

double* GPULodscores::gpu_init_lodscores() {
    double* tmp;
    
    CUDA_CALLANDTEST(cudaMalloc((void**)&tmp, sizeof(double) * (map->num_markers() - 1) * map->get_lodscore_count()));
    
    run_gpu_lodscoreinit_kernel(map->num_markers() - 1, map->get_lodscore_count(), tmp, map->get_lodscore_count());
    
    cudaThreadSynchronize();
    
    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess) {
        printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
        abort();
    }
    
    dev_lodscores = tmp;
    
    return tmp;
}

void GPULodscores::copy_to_gpu(DescentGraph& dg) {
    for(unsigned i = 0; i < unsigned(2 * ped->num_members() * map->num_markers()); ++i) {
        loc_state->dg->graph[i] = dg.get_bit(i);
    }
}

void GPULodscores::copy_from_gpu(DescentGraph& dg) {
    for(unsigned i = 0; i < unsigned(2 * ped->num_members() * map->num_markers()); ++i) {
        dg.set_bit(i, loc_state->dg->graph[i]);
    }
}

int GPULodscores::optimal_lodscore_threads() {
        
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
            
            run_gpu_lodscore_kernel(map->num_markers() - 1, threadcount[i], dev_state);
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

void GPULodscores::calculate(DescentGraph& dg) {

    block_until_finished();

    CUDA_CALLANDTEST(cudaMemcpy(dev_graph, dg.get_internal_ptr(), dg.get_internal_size(), cudaMemcpyHostToDevice));
    
    run_gpu_lodscore_kernel(map->num_markers() - 1, num_lodscore_threads, dev_state);
    
    ++count;
}

void GPULodscores::block_until_finished() {
    if(count != 0) {
        cudaThreadSynchronize();
                
        cudaError_t error = cudaGetLastError();
        if(error != cudaSuccess) {
            printf("CUDA kernel error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(error));
            abort();
        }
    }
}

void GPULodscores::get_results(LODscores* lod) {
    int num_lodscores = map->get_lodscore_count() * (map->num_markers() - 1);
    double* data = new double[num_lodscores];
    
    block_until_finished();
    
    CUDA_CALLANDTEST(cudaMemcpy(data, dev_lodscores, num_lodscores * sizeof(double), cudaMemcpyDeviceToHost));
    
    lod->set_trait_prob(trait_likelihood);
    lod->set_count(count);
    
    for(int i = 0; i < num_lodscores; ++i) {
        lod->set(i, data[i]);
    }
    
    free(data);
}

