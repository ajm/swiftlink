using namespace std;

#include <cstdio>
#include <iostream>
#include <fstream>

#include <omp.h>
#include <gsl/gsl_rng.h>

#include "random.h"

const gsl_rng_type* T;
gsl_rng** r;


void init_random() {    
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    //r = gsl_rng_alloc (T);
    
    r = new gsl_rng*[omp_get_max_threads()];
    
    for(int i = 0; i < omp_get_max_threads(); ++i) {
        r[i] = gsl_rng_alloc(T);
    }
}

void destroy_random() {
    //gsl_rng_free(r);
    
    for(int i = 0; i < omp_get_max_threads(); ++i) {
        gsl_rng_free(r[i]);
    }
    
    delete[] r;
}

unsigned int get_randomness() {
    unsigned int seed;
    fstream randfile;
    
    randfile.open ("/dev/urandom", ios::in);
    if(not randfile.is_open()) {
        fprintf(stderr, "error: could not open /dev/urandom\n");
        abort();
    }
      
    randfile.read((char*)&seed, sizeof(seed));
    
    randfile.close();
    
    return seed;
}

void seed_random_explicit(unsigned int seed) {
    //gsl_rng* ran = gsl_rng_alloc(T);
    //gsl_rng_set(ran, seed);
    
    //srandom(seed);
    
    for(int i = 0; i < omp_get_max_threads(); ++i) {
        //gsl_rng_set(r[i], random());
        //gsl_rng_set(r[i], gsl_rng_get(ran));
        //gsl_rng_set(r[i], get_randomness());
        
        unsigned int s = get_randomness();
        gsl_rng_set(r[i], s);
        
        printf("seed %d = %u\n", i, s);
    }
    
    printf("generator type: %s\n", gsl_rng_name(r[0]));
    //printf("seed = %u\n", seed);
    
    //gsl_rng_free(ran);
}

void seed_random_implicit() {
    seed_random_explicit(get_randomness());
}

double get_random() {
    return gsl_rng_uniform(r[omp_get_thread_num()]);
}

int get_random_int(int limit) {
    return gsl_rng_uniform_int(r[omp_get_thread_num()], limit);
}

