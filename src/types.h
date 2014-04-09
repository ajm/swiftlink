#ifndef LKG_TYPES_H_
#define LKG_TYPES_H_

#include <cstdlib>

#include "genotype.h"
#include "trait.h"
#include "defaults.h"


const unsigned int UNKNOWN_PARENT = ~0u;
const unsigned int UNKNOWN_ID = UNKNOWN_PARENT;
const unsigned int DEBUG_FP_PRECISION = 3;
const double DBL_RAND_MAX = static_cast<double>(RAND_MAX);

enum parentage { 
    MATERNAL,
    PATERNAL
};

enum sex {
    UNSEXED,
    MALE,
    FEMALE
};

enum affection {
    UNKNOWN_AFFECTION,
    UNAFFECTED,
    AFFECTED
};

enum simple_disease_model {
    AUTOSOMAL_RECESSIVE,
    AUTOSOMAL_DOMINANT
};

// used in a few places following nomenclature of
// Corman, Leiserson, Rivest & Stein 2nd Ed. depth-first search
enum {
    WHITE,
    GREY,
    BLACK
};

typedef int meiosis_indicator_t;
typedef enum parentage allele_t;

string gender_str(enum sex s);
string affection_str(enum affection a);
string parent_str(enum parentage p);

struct mcmc_options {
    bool verbose;

    int burnin;
    int iterations;
    int si_iterations;
    int scoring_period;
    int lodscores;
    int peelopt_iterations;
    
    double lsampler_prob;
    
    // parallelism
    int thread_count;
    bool use_gpu;
    
    // things precalculated or stored in files
    string peelseq_filename;
    string random_filename;

    // elod options
    double elod_frequency;
    double elod_penetrance[3];
    double elod_marker_separation;

    // mc3 options
    int mc3_number_of_chains;
    int mc3_exchange_period;
    double mc3_temperature;

    mcmc_options() :
        verbose(false),
        burnin(DEFAULT_BURNIN_ITERATIONS),
        iterations(DEFAULT_MCMC_ITERATIONS),
        si_iterations(DEFAULT_SEQUENTIALIMPUTATION_RUNS),
        scoring_period(DEFAULT_SCORING_PERIOD),
        lodscores(DEFAULT_LODSCORES),
        peelopt_iterations(DEFAULT_PEELOPT_ITERATIONS),
        lsampler_prob(DEFAULT_LSAMPLER_PROB),
        thread_count(DEFAULT_THREAD_COUNT),
        use_gpu(false),
        peelseq_filename(""),
        random_filename(""),
        elod_frequency(0.0),
        elod_marker_separation(0.0), 
        mc3_number_of_chains(0), 
        mc3_exchange_period(0),
        mc3_temperature(0.0) {
        
            elod_penetrance[0] = elod_penetrance[1] = elod_penetrance[2] = 0.0;
        }
};


#endif

