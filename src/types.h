#ifndef LKG_TYPES_H_
#define LKG_TYPES_H_

#include <cstdlib>
#include <vector>
#include <sstream>

#include "trait.h"
#include "defaults.h"


const unsigned int UNKNOWN_PARENT = ~0u;
const unsigned int UNKNOWN_ID = UNKNOWN_PARENT;
const unsigned int DEBUG_FP_PRECISION = 3;
const double DBL_RAND_MAX = static_cast<double>(RAND_MAX);

enum parentage { 
    MATERNAL,
    PATERNAL,
    NONE
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

    // mcmc
    int burnin;
    int iterations;
    int si_iterations;
    int scoring_period;
    int mcmc_runs;

    // linkage
    int lodscores;
    int peelopt_iterations;
    
    double lsampler_prob;
    
    // parallelism
    int thread_count;
    bool use_gpu;
    
    // things precalculated or stored in files
    string peelseq_filename;
    string random_filename;
    string exchange_filename;

    bool affected_only;
    bool sex_linked;

    // elod options
    bool elod;
    double elod_frequency;
    vector<double> elod_penetrance;
    double elod_marker_separation;
    int elod_replicates;

    // mc3 options
    int mc3;
    int mc3_number_of_chains;
    int mc3_exchange_period;
    vector<double> mc3_temperatures;

    mcmc_options() :
        verbose(false),
        burnin(DEFAULT_MCMC_BURNIN),
        iterations(DEFAULT_MCMC_ITERATIONS),
        si_iterations(DEFAULT_SEQUENTIALIMPUTATION_RUNS),
        scoring_period(DEFAULT_MCMC_SCORING_PERIOD),
        mcmc_runs(DEFAULT_MCMC_RUNS),
        lodscores(DEFAULT_LODSCORES),
        peelopt_iterations(DEFAULT_PEELOPT_ITERATIONS),
        lsampler_prob(DEFAULT_LSAMPLER_PROB),
        thread_count(DEFAULT_THREAD_COUNT),
        use_gpu(false),
        peelseq_filename(""),
        random_filename(""),
        exchange_filename(""),
        affected_only(false),
        sex_linked(false),
        elod(false),
        elod_frequency(DEFAULT_ELOD_FREQUENCY),
        elod_penetrance(DEFAULT_ELOD_PENETRANCE, DEFAULT_ELOD_PENETRANCE + 3),
        elod_marker_separation(DEFAULT_ELOD_SEPARATION),
        elod_replicates(DEFAULT_ELOD_REPLICATES),
        mc3(false),
        mc3_number_of_chains(DEFAULT_MCMC_CHAINS), 
        mc3_exchange_period(DEFAULT_MCMC_EXCHANGE_PERIOD),
        mc3_temperatures() {}

    string debug_string() {
        stringstream ss;

        ss << "elod_penetrance: ";
        for(int i = 0; i < int(elod_penetrance.size()); ++i)
            ss << elod_penetrance[i] << " ";
        ss << "\n";

        ss << "mc3_temperatures: ";
        for(int i = 0; i < int(mc3_temperatures.size()); ++i)
            ss << mc3_temperatures[i] << " ";
        ss << "\n";

        return ss.str();
    }
};

#endif

