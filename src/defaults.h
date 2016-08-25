#ifndef LKG_DEFAULTS_H_
#define LKG_DEFAULTS_H_

const double DEFAULT_ELOD_FREQUENCY         = 0.00001;
const double DEFAULT_ELOD_PENETRANCE[3]     = { 0.0, 0.0, 1.0 };
const double DEFAULT_ELOD_SEPARATION        = 0.05;
const int DEFAULT_ELOD_REPLICATES           = 1000000;

const int DEFAULT_MCMC_CHAINS               = 1;
const int DEFAULT_MCMC_ITERATIONS           = 90000;
const int DEFAULT_MCMC_BURNIN               = 10000;
const int DEFAULT_MCMC_EXCHANGE_PERIOD      = 10;
const int DEFAULT_MCMC_SCORING_PERIOD       = 10;
const int DEFAULT_MCMC_RUNS                 = 1;

#define DEFAULT_RESULTS_FILENAME "swiftlink.out"
#define DEFAULT_CODA_PREFIX "coda"

const int DEFAULT_SEQUENTIALIMPUTATION_RUNS = 1000;

const int DEFAULT_THREAD_COUNT              = 1;
const int DEFAULT_LODSCORES                 = 5;
const int DEFAULT_PEELOPT_ITERATIONS        = 1000000;
const double DEFAULT_LSAMPLER_PROB          = 0.5;

#endif

