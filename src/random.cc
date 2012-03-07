using namespace std;

#include <cstdio>
#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>

#include "random.h"

const gsl_rng_type* T;
gsl_rng* r;


void init_random() {    
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);    
}

void destroy_random() {
    gsl_rng_free(r);
}

void seed_random_explicit(unsigned int seed) {
    gsl_rng_set(r, seed);

    printf ("generator type: %s\n", gsl_rng_name(r));
    printf ("seed = %u\n", seed);
}

void seed_random_implicit() {
    unsigned int seed;
    fstream randfile;
    
    randfile.open ("/dev/urandom", ios::in);
    if(not randfile.is_open()) {
        fprintf(stderr, "error: could not open /dev/urandom\n");
        abort();
    }
      
    randfile.read((char*)&seed, sizeof(seed));
    
    randfile.close();
    
    
    seed_random_explicit(seed);
}

double get_random() {
    return gsl_rng_uniform(r);
}

int get_random_int(int limit) {
    return gsl_rng_uniform_int(r, limit);
}

