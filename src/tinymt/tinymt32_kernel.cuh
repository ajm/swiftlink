#ifndef TINYMT32_KERNEL_CUH
#define TINYMT32_KERNEL_CUH
/**
 * @file tinymt32_kernel.cuh
 *
 * @brief CUDA implementation of TinyMT32.
 *
 * This is CUDA implementation of TinyMT32 pseudo-random number generator.
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (The University of Tokyo)
 *
 * Copyright (C) 2011 Mutsuo Saito, Makoto Matsumoto,
 * Hiroshima University and The University of Tokyo.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see LICENSE.txt
 */
#include <cuda.h>
#include <stdint.h>

// i moved the tinymt32_status_t definition to here
// so I could use it in C++  --ajm
#include "tinymt32_host.h"

/* =====================================
   DEFINITIONS FOR USERS
   ===================================== */
/**
 * TinyMT structure
 * mat1, mat2, tmat must be set from tinymt32dc output before initialize
 */
/*
typedef struct TINYMT32_STATUS_T {
    uint32_t status[4];
    uint32_t mat1;
    uint32_t mat2;
    uint32_t tmat;
} tinymt32_status_t;
*/

/**
 * Initialize TinyMT structure by seed and parameters.
 *
 * This function must be called before tinymt32_uint32(),
 * tinymt32_single(), tinymt32_single12().
 * mat1, mat2, tmat in tinymt must be set before this call.
 *
 * @param tinymt TinyMT structure to be initialized.
 * @param seed seed of randomness.
 */
__device__ void tinymt32_init(tinymt32_status_t* tinymt, uint32_t seed);

/**
 * Generate 32bit unsigned integer r (0 <= r < 2<sup>32</sup>)
 *
 * @param tinymt TinyMT structure
 * @return 32bit unsigned integer
 */
__device__ uint32_t tinymt32_uint32(tinymt32_status_t * tinymt);

/**
 * Generate single precision floating point number r (0.0 <= r < 1.0)
 *
 * @param tinymt TinyMT structure
 * @return single precision floating point number
 */
__device__ float tinymt32_single(tinymt32_status_t * tinymt);

/**
 * Generate single precision floating point number r (1.0 <= r < 2.0).
 *
 * This function may little bit faster than tinymt32_single().
 *
 * @param tinymt TinyMT structure
 * @return single precision floating point number
 */
__device__ float tinymt32_single12(tinymt32_status_t * tinymt);

/* =====================================
   DEFINITIONS FOR INTERNAL USE
   ===================================== */
#define TINYMT32_SHIFT0 1
#define TINYMT32_SHIFT1 10
#define TINYMT32_MIN_LOOP 8
#define TINYMT32_PRE_LOOP 8
#define TINYMT32_MASK 0x7fffffffU
#define TINYMT32_SINGLE_MASK 0x3f800000U

__device__ void tinymt32_next(tinymt32_status_t * tinymt);
__device__ uint32_t tinymt32_temper(tinymt32_status_t * tinymt);

/* =====================================
   FUNCTIONS
   ===================================== */
/**
 * The state transition function.
 * @param tinymt the internal state.
 */
__device__ void tinymt32_next(tinymt32_status_t * tinymt)
{
    uint32_t y = tinymt->status[3];
    uint32_t x = (tinymt->status[0] & TINYMT32_MASK)
	^ tinymt->status[1] ^ tinymt->status[2];
    x ^= (x << TINYMT32_SHIFT0);
    y ^= (y >> TINYMT32_SHIFT0) ^ x;
    tinymt->status[0] = tinymt->status[1];
    tinymt->status[1] = tinymt->status[2];
    tinymt->status[2] = x ^ (y << TINYMT32_SHIFT1);
    tinymt->status[3] = y;
    if (y & 1) {
	tinymt->status[1] ^= tinymt->mat1;
	tinymt->status[2] ^= tinymt->mat2;
    }
}

/**
 * The tempering function.
 *
 * This function improves the equidistribution property of
 * outputs.
 * @param tinymt the internal state.
 * @return tempered output
 */
__device__ uint32_t tinymt32_temper(tinymt32_status_t * tinymt)
{
    uint32_t t0, t1;
    t0 = tinymt->status[3];
    t1 = tinymt->status[0] + (tinymt->status[2] >> 8);
    t0 ^= t1;
    if (t1 & 1) {
	t0 ^= tinymt->tmat;
    }
    return t0;
}

__device__ uint32_t tinymt32_uint32(tinymt32_status_t * tinymt)
{
    tinymt32_next(tinymt);
    return tinymt32_temper(tinymt);
}

__device__ float tinymt32_single12(tinymt32_status_t * tinymt)
{
    uint32_t t0;
    tinymt32_next(tinymt);
    t0 = tinymt32_temper(tinymt);
    t0 = t0 >> 9;
    t0 ^= TINYMT32_SINGLE_MASK;
    return __int_as_float(t0);
}

__device__ float tinymt32_single(tinymt32_status_t * tinymt)
{
    return tinymt32_single12(tinymt) - 1.0f;
}

__device__ void tinymt32_init(tinymt32_status_t* tinymt, uint32_t seed) {
    tinymt->status[0] = seed;
    tinymt->status[1] = tinymt->mat1;
    tinymt->status[2] = tinymt->mat2;
    tinymt->status[3] = tinymt->tmat;
    for (int i = 1; i < TINYMT32_MIN_LOOP; i++) {
	tinymt->status[i & 3] ^= i + 1812433253U *
	    (tinymt->status[(i - 1) & 3]
	     ^ (tinymt->status[(i - 1) & 3] >> 30));
    }
    if ((tinymt->status[0] & TINYMT32_MASK) == 0 &&
	tinymt->status[1] == 0 &&
	tinymt->status[2] == 0 &&
	tinymt->status[3] == 0) {
	tinymt->status[0] = 'T';
	tinymt->status[1] = 'I';
	tinymt->status[2] = 'N';
	tinymt->status[3] = 'Y';
    }
    for (int i = 0; i < TINYMT32_PRE_LOOP; i++) {
	tinymt32_next(tinymt);
    }
}

#undef TINYMT32_SHIFT0
#undef TINYMT32_SHIFT1
#undef TINYMT32_MIN_LOOP
#undef TINYMT32_PRE_LOOP
#undef TINYMT32_MASK
#undef TINYMT32_SINGLE_MASK

#endif

