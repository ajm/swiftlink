#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#if __GNUC__
#if __x86_64__
#error "64 bit!"
#else
//#error "32 bit!"
//#include "mt19937-32/mt19937ar.c"
#include "dc.h"
#endif
#endif

//#include "types.h"
#include "random.h"



static mt_struct **mtss;
static int mtss_count;


void seed_random(unsigned int seed) {
    //srandom(seed);
    //init_genrand(seed);
    int i;
    
    srandom(seed);
    
    // if it has been initialised before
    if(mtss) {
        free_mt_struct_array(mtss, mtss_count);
    }
    
    mtss_count = omp_get_max_threads() - 1;
    printf("omp_get_num_threads() = %d\n", mtss_count + 1);
    
    // taken from new_example3.c, dcmt 0.6.1
    mtss = get_mt_parameters_st(32, 521, 0, mtss_count, 4172, &mtss_count);
    
    if(!mtss) {
        fprintf(stderr, "error: failed to initialise parallel prng (%s:%d)\n", __FILE__, __LINE__);
        abort();
    }
    
    printf("generated %d prng\n", mtss_count);
    
    // seed
    for(i = 0; i < mtss_count; ++i) {
        sgenrand_mt(random(), mtss[i]);
    }
}

double get_random() {
    //return random() / DBL_RAND_MAX;
    //return genrand_real1();
    return genrand_mt(mtss[omp_get_thread_num()]) * (1.0/4294967295.0);
}

int get_random_int(int limit) {
    //return get_random() % limit;
    //return genrand_int32() % limit;
    return genrand_mt(mtss[omp_get_thread_num()]) % limit;
}

void destroy_random() {
    free_mt_struct_array(mtss, mtss_count);
    mtss_count = 0;
}

