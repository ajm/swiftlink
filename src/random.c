#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <float.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

/*
#if __GNUC__
#if __x86_64__
#error "64 bit!"
#else
//#error "32 bit!"
//#include "mt19937-32/mt19937ar.c"
#include "dc.h"
#endif
#endif
*/

//#include "dc.h"
//#include "mt64.h"
#include "mt19937-32/mt19937ar.c"


//#include "types.h"
#include "random.h"

#include <gsl/gsl_rng.h>

const gsl_rng_type* T;
gsl_rng* r;

/*
static mt_struct **mtss;
static int mtss_count;

static int fd;
*/
/*
unsigned int get_randomness() {
    unsigned int seed;
    
    if(read(fd, &seed, sizeof(seed)) != sizeof(seed)) {
        abort();
    }
    
    return seed;
}
*/
void seed_random(unsigned int seed) {
    
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    gsl_rng_set(r, seed);

    printf ("generator type: %s\n", gsl_rng_name (r));
    printf ("seed = %lu\n", gsl_rng_default_seed);

    return;
    
    
    //srandom(seed);
    
    // XXX
    //init_genrand64(seed);
    
    //seed = 1206534962;
    //seed = 327423963;
    //seed = 2316696976;
    
    printf("seed = %u\n", seed);
    init_genrand(seed); //5489);
    
    
    /*
    int i;
    
    
    
    fd = open("/dev/urandom", O_RDONLY);
    
    
    
    srandom(seed);
    
    printf("seed = %u\n", seed);
    
    
    // if it has been initialised before
    if(mtss) {
        free_mt_struct_array(mtss, mtss_count);
    }
    
    mtss_count = omp_get_max_threads() - 1;
    printf("omp_get_num_threads() = %d\n", mtss_count + 1);
    
    // taken from new_example3.c, dcmt 0.6.1
    //mtss = get_mt_parameters_st(32, 521, 0, mtss_count, 4172, &mtss_count);
    mtss = get_mt_parameters_st(32, 1279, 0, mtss_count, get_randomness(), &mtss_count);
    
    if(!mtss) {
        fprintf(stderr, "error: failed to initialise parallel prng (%s:%d)\n", __FILE__, __LINE__);
        abort();
    }
    
    printf("generated %d prng\n", mtss_count);
    
    // seed
    for(i = 0; i < mtss_count; ++i) {
        //sgenrand_mt(random(), mtss[i]);
        sgenrand_mt(get_randomness(), mtss[i]);
    }
    
    close(fd);
    */
}

double get_random() {

    return gsl_rng_uniform(r);

    //return (random() * (1.0/4294967296.0)) + (1.0/4294967296.0);
    
    // XXX
    //return ((genrand64_int64() >> 11) * (1.0/9007199254740992.0)) + (1.0/9007199254740992.0);
    
    return genrand_real3();
    
    //double r = (double) genrand_mt(mtss[omp_get_thread_num()]);
    //return r * (1.0/4294967295.0);
    
    // generate random number in the interval [0,1]
    //double r = genrand_mt(mtss[omp_get_thread_num()]) * (1.0/4294967295.0);
    
    /*
    // generate random number in the interval [0,1)
    double r = (genrand_mt(mtss[omp_get_thread_num()]) * (1.0/4294967296.0)) + (1.0/4294967296.0);
    
    if((r <= 0.0) || (r > 1.0)) {
        fprintf(stderr, "%f\n", r);
        abort();
    }
    
    return r;
    */
}

unsigned int get_random_raw() {
    return genrand_int32();
}

int get_random_int(int limit) {

    return gsl_rng_uniform_int(r, limit);

    //return random() % limit;
    // XXX
    //return genrand64_int64() % limit;
    
    return genrand_int32() % limit;
    //return genrand_mt(mtss[omp_get_thread_num()]) % limit;
}

void destroy_random() {
    //free_mt_struct_array(mtss, mtss_count);
    //mtss_count = 0;
}

