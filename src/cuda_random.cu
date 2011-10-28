#include "cuda_common.h"
#include "tinymt/tinymt32_kernel.cuh"

//#define RANDOM_USE_MT 1


__device__ float get_curand_random(struct gpu_state* state) {
    //return random() / (float) RAND_MAX;
    curandState localstate;
    //int offset;
    float tmp;
    
    //offset = (blockIdx.x * 256) + threadIdx.x;
    
    localstate = state->randstates[blockIdx.x];
    
    tmp = curand_uniform(&localstate);
    
    state->randstates[blockIdx.x] = localstate;

    return tmp;
}

__device__ float get_mt_random(struct gpu_state* state) {
    
    return tinymt32_single(&state->randmt[blockIdx.x]);
}

__device__ float get_random(struct gpu_state* state) {

    //return 0.5;
#ifdef RANDOM_USE_MT
    return get_mt_random(state);
#else
    return get_curand_random(state);
#endif
}

__global__ void init_curand_kernel(curandState* states, long int* seeds) {
    int id = blockIdx.x; // (blockIdx.x * 256) + threadIdx.x;
    
    //printf("%d %ld\n", id, seeds[id]);
    
    curand_init(seeds[id], 0, 0, &states[id]);
    //curand_init(id, 0, 0, &states[id]);
    //curand_init(1234, id, 0, &states[id]);
}

__global__ void init_mt_kernel(tinymt32_status_t* states, uint32_t* params, uint32_t* seeds) {
    int id = blockIdx.x; //(blockIdx.x * 256) + threadIdx.x;
    
    states->mat1 = params[id * 3];
    states->mat2 = params[(id * 3) + 1];
    states->tmat = params[(id * 3) + 2];

    tinymt32_init(&(states[id]), seeds[id]);
}

void run_gpu_curand_init_kernel(int numblocks, curandState* states, long int* seeds) {
    init_curand_kernel<<<numblocks, 1>>>(states, seeds);
}

void run_gpu_tinymt_init_kernel(int numblocks, tinymt32_status_t* states, uint32_t* params, uint32_t* seeds) {
    init_mt_kernel<<<numblocks, 1>>>(states, params, seeds);
}

