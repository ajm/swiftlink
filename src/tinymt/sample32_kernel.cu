/**
 * @file sample32_kernel.cu
 *
 * @brief Sample Program for CUDA 4.0
 *
 * sample kernel program which call tinymt32_single().
 */
#include <stdint.h>
#include "tinymt32_kernel.cuh"

/**
 * kernel function.
 * This function generates single precision floating point numbers
 * and sum up the numbers for each thread.
 *
 * -# set parameters
 * -# initialize the structure
 * -# generate
 *
 * @param param_array parameters array of TinyMT
 * @param sum_array the array of the result
 * @param size number of output data sum up for each thread.
 */
__global__ void sample_single_sum_kernel(uint32_t * param_array,
					 float* sum_array,
					 int size)
{
    const int tid = blockDim.x * blockIdx.x + threadIdx.x;
    const uint32_t *p = &param_array[tid * 3];
    tinymt32_status_t tinymt32;
    float sum = 0.0f;

    // set parameters
    tinymt32.mat1 = *p;
    tinymt32.mat2 = *(p + 1);
    tinymt32.tmat = *(p + 2);

    // initialize
    tinymt32_init(&tinymt32, tid);

    for (int i = 0; i < size; i++) {
	// generate
	sum += tinymt32_single(&tinymt32);
	//sum += tinymt32_single12(&tinymt32);
	//sum += tinymt32_uint32(&tinymt32);
    }
    sum_array[tid] = sum;
}

