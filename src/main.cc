/*
 * the only thing this file does is handle command line arguments
 * and then call another main function
 */
 
#define _GNU_SOURCE 1
#include <fenv.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <climits>
#include <unistd.h>
#include <getopt.h>

#include "types.h"
#include "defaults.h"
#include "linkage_program.h"
#include "elod.h"
#include "omp_facade.h"

//#include "test_program.h"
//#include "haplotype_program.h"

/*
enum analysistype {
    LINKAGE,
    HAPLOTYPE,
    TESTING
};
*/

char* mapfile = NULL;
char* pedfile = NULL;
char* datfile = NULL;
char* outfile = DEFAULT_RESULTS_FILENAME;

struct mcmc_options options;

//enum analysistype analysis = LINKAGE;


void _usage(char *prog) {
	fprintf(stderr,
"Usage: %s [OPTIONS] -p pedfile -m mapfile -d datfile\n"
"       %s [OPTIONS] -p pedfile --elod\n"
"\n"
"Input files:\n"
"  -p pedfile, --pedigree=pedfile\n"
"  -m mapfile, --map=mapfile\n"
"  -d datfile, --dat=datfile\n"
"\n"
"Output files:\n"
"  -o outfile, --output=outfile            (default = '%s')\n"
"\n"
"MCMC options:\n"
"  -i NUM,     --iterations=NUM            (default = %d)\n"
"  -b NUM,     --burnin=NUM                (default = %d)\n"
"  -s NUM,     --sequentialimputation=NUM  (default = %d)\n"
"  -x NUM,     --scoringperiod=NUM         (default = %d)\n"
"  -l FLOAT,   --lsamplerprobability=FLOAT (default = %.1f)\n"
"  -n NUM,     --lodscores=NUM             (default = %d)\n"
"  -R NUM,     --runs=NUM                  (default = %d)\n"
"\n"
"MCMC diagnostic options:\n"
"  -C,         --coda\n"
"  -P PREFIX,  --codaprefix=PREFIX         (default = %s)\n"
"\n"
//"Metropolis-coupled MCMC options:\n"
//"  -M,         --mcmcmc\n"
//"  -z NUM,     --chains=NUM                (default = %d)\n"
//"  -y NUM,     --exchangeperiod=NUM        (default = %d)\n"
//"  -t FLOAT,FLOAT,... --temperatures=FLOAT,FLOAT,...\n"
//"\n"
"ELOD options:\n"
"  -e          --elod\n"
"  -f FLOAT    --frequency=FLOAT           (default = %.4f)\n"
"  -w FLOAT    --separation=FLOAT          (default = %.4f)\n"
"  -k FLOAT,FLOAT,FLOAT --penetrance=FLOAT,FLOAT,FLOAT(default = %.2f,%.2f,%.2f)\n"
"  -u NUM      --replicates=NUM            (default = %d)\n"
"\n"
"Runtime options:\n"
"  -c NUM,     --cores=NUM                 (default = %d)\n"
#ifndef USE_CUDA
"  -g,         --gpu                       [UNAVAILABLE, COMPILED WITHOUT CUDA]\n"
#else
"  -g,         --gpu\n"
#endif
"\n"
"Misc:\n"
"  -X,         --sexlinked\n"
"  -a,         --affectedonly\n"
"  -q NUM,     --peelseqiter=NUM           (default = %d)\n"
"  -r seedfile,--randomseeds=seedfile\n"
"  -v,         --verbose\n"
"  -h,         --help\n"
"\n", 
prog, 
prog,
DEFAULT_RESULTS_FILENAME, 
DEFAULT_MCMC_ITERATIONS,
DEFAULT_MCMC_BURNIN,
DEFAULT_SEQUENTIALIMPUTATION_RUNS,
DEFAULT_MCMC_SCORING_PERIOD,
DEFAULT_LSAMPLER_PROB,
DEFAULT_LODSCORES,
DEFAULT_MCMC_RUNS,
DEFAULT_CODA_PREFIX,
//DEFAULT_MCMC_CHAINS,
//DEFAULT_MCMC_EXCHANGE_PERIOD,
DEFAULT_ELOD_FREQUENCY,
DEFAULT_ELOD_SEPARATION,
DEFAULT_ELOD_PENETRANCE[0],
DEFAULT_ELOD_PENETRANCE[1],
DEFAULT_ELOD_PENETRANCE[2],
DEFAULT_ELOD_REPLICATES,
DEFAULT_THREAD_COUNT,
DEFAULT_PEELOPT_ITERATIONS
);
}

bool str2int(int &i, char* s) {
    char* end;
    long l;
    errno = 0;
    l = strtol(s, &end, 10);
    
    if (((errno == ERANGE) and (l == LONG_MAX)) or (l > INT_MAX)) {
        return false;
    }
    if (((errno == ERANGE) and (l == LONG_MIN)) or (l < INT_MIN)) {
        return false;
    }
    if ((*s == '\0') or (*end != '\0')) {
        return false;
    }
    
    i = l;
    return true;
}

bool str2float(double &i, char *s) {
    char* end;
    double d;
    errno = 0;
    
    d = strtod(s, &end);
    
    // overflow
    if((errno == ERANGE) and (d == HUGE_VAL)) {
        return false;
    }
    
    // underflow
    if((errno == ERANGE) and (d == 0.0)) {
        return false;
    }
        
    if ((*s == '\0') or (*end != '\0')) {
        return false;
    }
    
    i = d;
    return true;
}

bool csv2vec(vector<double>& v, char *s) {
    char* tmp;
    double d;

    v.clear();

    while((tmp = strsep(&s, ",")) != NULL) {
        if(not str2float(d, tmp)) {
            return false;
        }

        v.push_back(d);
    }

    return true;
}

void _handle_args(int argc, char **argv) {
	extern char *optarg;
    extern int optopt;
	int ch;
	
	static struct option long_options[] = 
	    {
            {"affectedonly",        no_argument,        0,      'a'},
            {"burnin",              required_argument,  0,      'b'},
            {"cores",               required_argument,  0,      'c'},
            {"dat",                 required_argument,  0,      'd'},
            {"elod",                no_argument,        0,      'e'},
            {"frequency",           required_argument,  0,      'f'},
            {"gpu",                 no_argument,        0,      'g'},
            {"help",                no_argument,        0,      'h'},
            {"iterations",          required_argument,  0,      'i'},
            //{"exchangefile",        required_argument,  0,      'j'},
            {"penetrance",          required_argument,  0,      'k'},
            {"lsamplerprobability", required_argument,  0,      'l'},
            {"map",                 required_argument,  0,      'm'},
            {"lodscores",           required_argument,  0,      'n'},
            {"output",              required_argument,  0,      'o'},
            {"pedigree",            required_argument,  0,      'p'},
            {"peelseqiter",         required_argument,  0,      'q'},
            {"randomseeds",         required_argument,  0,      'r'},
            {"sequentialimputation",required_argument,  0,      's'},
            //{"temperatures",        required_argument,  0,      't'},
            {"replicates",          required_argument,  0,      'u'},
	        {"verbose",             no_argument,        0,      'v'},
            {"separation",          required_argument,  0,      'w'},
            {"scoringperiod",       required_argument,  0,      'x'},
            //{"exchangeperiod",      required_argument,  0,      'y'},
            //{"chains",              required_argument,  0,      'z'},
            //{"mcmcmc",              no_argument,        0,      'M'},
            {"sexlinked",           no_argument,        0,      'X'},
            {"runs",                required_argument,  0,      'R'},
            {"coda",                no_argument,        0,      'C'},
            {"codaprefix",          required_argument,  0,      'P'},
            {0, 0, 0, 0}
	    };
    
    int option_index = 0;
    
	while ((ch = getopt_long(argc, argv, 
                    //":p:d:m:o:i:b:s:l:c:x:q:r:n:vhcgz:y:t:ew:k:f:u:j:aMX", 
                    ":p:d:m:o:i:b:s:l:c:x:q:r:n:vhcgew:k:f:u:aXR:CP:",
                    long_options, &option_index)) != -1) {
		switch (ch) {
			case 'p':
				pedfile = optarg;
				break;
				
			case 'm':
				mapfile = optarg;
				break;
				
			case 'd':
				datfile = optarg;
				break;
				
			case 'v':
				options.verbose = true;
				break;
				
            case 'o':
                outfile = optarg;
                break;
                
            case 'i':
                if(not str2int(options.iterations, optarg)) {
                    fprintf(stderr, "%s: option '-i' requires an int as an argument ('%s' given)\n", argv[0], optarg);
                    exit(EXIT_FAILURE);
                }
                if(options.iterations <= 0) {
                    fprintf(stderr, "%s: mcmc iterations must be positive (%d given)\n", argv[0], options.iterations);
                    exit(EXIT_FAILURE);
                }
                break;
                
            case 'b':
                if(not str2int(options.burnin, optarg)) {
                    fprintf(stderr, "%s: option '-b' requires an int as an argument ('%s' given)\n", argv[0], optarg);
                    exit(EXIT_FAILURE);
                }
                if(options.burnin < 0) {
                    fprintf(stderr, "%s: mcmc burnin must be positive (%d given)\n", argv[0], options.burnin);
                    exit(EXIT_FAILURE);
                }
                break;
                
            case 's':
                if(not str2int(options.si_iterations, optarg)) {
                    fprintf(stderr, "%s: option '-s' requires an int as an argument ('%s' given)\n", argv[0], optarg);
                    exit(EXIT_FAILURE);
                }
                if(options.si_iterations < 1) {
                    fprintf(stderr, "%s: number of sequential imputation runs must be greater than zero (%d given)\n", argv[0], options.si_iterations);
                    exit(EXIT_FAILURE);
                }
                break;
                
            case 'x':
                if(not str2int(options.scoring_period, optarg)) {
                    fprintf(stderr, "%s: option '-x' requires an int as an argument ('%s' given)\n", argv[0], optarg);
                    exit(EXIT_FAILURE);
                }
                if(options.scoring_period <= 0) {
                    fprintf(stderr, "%s: number of sequential imputation runs must be greater than 0 (%d given)\n", argv[0], options.scoring_period);
                    exit(EXIT_FAILURE);
                }
                break;
                
            case 'c':
                if(not str2int(options.thread_count, optarg)) {
                    fprintf(stderr, "%s: option '-c' requires an int as an argument ('%s' given)\n", argv[0], optarg);
                    exit(EXIT_FAILURE);
                }
                if(options.thread_count < 1) {
                    fprintf(stderr, "%s: thread count must be at least one (%d given)\n", argv[0], options.thread_count);
                    exit(EXIT_FAILURE);
                }
                break;
                
            case 'l':
                if(not str2float(options.lsampler_prob, optarg)) {
                    fprintf(stderr, "%s: option '-l' requires a float (>=0, <= 1) as an argument ('%s' given)\n", argv[0], optarg);
                    exit(EXIT_FAILURE);
                }
                if((options.lsampler_prob < 0.0) or (options.lsampler_prob > 1.0)) {
                    fprintf(stderr, "%s: option '-l' requires an argument >= 0 and <= 1 (%.3f given)\n", argv[0], options.lsampler_prob);
                    exit(EXIT_FAILURE);
                }
                break;
                
            case 'n':
                if(not str2int(options.lodscores, optarg)) {
                    fprintf(stderr, "%s: option '-n' requires an int as an argument ('%s' given)\n", argv[0], optarg);
                    exit(EXIT_FAILURE);
                }
                if(options.lodscores <= 0) {
                    fprintf(stderr, "%s: lodscores must be positive (%d given)\n", argv[0], options.lodscores);
                    exit(EXIT_FAILURE);
                }
                break;
                
            case 'g':
#ifdef USE_CUDA
                options.use_gpu = true;
                break;
#else
                fprintf(stderr, "Error: SwiftLink was compiled without CUDA support, exiting...\n");
                exit(EXIT_FAILURE);
#endif
                
            case 'q':
                if(not str2int(options.peelopt_iterations, optarg)) {
                    fprintf(stderr, "%s: option '-q' requires an int as an argument ('%s' given)\n", argv[0], optarg);
                    exit(EXIT_FAILURE);
                }
                if(options.peelopt_iterations <= 0) {
                    fprintf(stderr, "%s: peeling sequence iterations must be positive (%d given)\n", argv[0], options.lodscores);
                    exit(EXIT_FAILURE);
                }
                // if this is too small, then you can get floating point exceptions
                // due to the summation of the sizes of the matrices being too big
                if(options.peelopt_iterations < 10000) {
                    options.peelopt_iterations = 10000;
                }
                break;
                
            case 'r':
                options.random_filename = string(optarg);
                break;
            
            /*
		    case 'a':
		        if(strcmp(optarg, "linkage") == 0) {
		            analysis = LINKAGE;
		        }
		        else if(strcmp(optarg, "haplotype") == 0) {
		            analysis = HAPLOTYPE;
		        }
		        else if(strcmp(optarg, "test") == 0) {
		            analysis = TESTING;
		        }
		        else {
		            fprintf(stderr, "%s: option '-a' can only accept 'linkage' or "
                                    "'haplotype' as arguments ('%s' given)\n",
                                    argv[0], optarg);
		            exit(EXIT_FAILURE);
		        }
		        break;
		    */
			case 'h':
				_usage(argv[0]);
				exit(EXIT_SUCCESS);
	        /*
            case 'j':
                options.exchange_filename = string(optarg);
                break;

            case 'z':
                if(not str2int(options.mc3_number_of_chains, optarg)) {
                    fprintf(stderr, "%s: option '-z' requires an int as an argument ('%s' given)\n", argv[0], optarg);
                    exit(EXIT_FAILURE);
                }
                if(options.mc3_number_of_chains < 1) {
                    fprintf(stderr, "%s: number of Markov chains must be greater than zero ('%d' given)\n", argv[0], options.mc3_number_of_chains);
                    exit(EXIT_FAILURE);
                }
                break;

            case 'y':
                if(not str2int(options.mc3_exchange_period, optarg)) {
                    fprintf(stderr, "%s: option '-y' requires an int as an argument ('%s' given)\n", argv[0], optarg);
                    exit(EXIT_FAILURE);
                }
                if(options.mc3_exchange_period < 1) {
                    fprintf(stderr, "%s: exchange period must be greater than zero ('%d' given)\n", argv[0], options.mc3_exchange_period);
                    exit(EXIT_FAILURE);
                }
                break;

            case 't':
                if(not csv2vec(options.mc3_temperatures, optarg)) {
                    fprintf(stderr, "%s: temperatures must be comma delimited floats from 0.0 - 1.0 inclusive, e.g.: 1.0,0.9,0.8,0.7 ('%s' given)\n", argv[0], optarg);
                    exit(EXIT_FAILURE);
                }
                
                for(unsigned i = 0; i < options.mc3_temperatures.size(); ++i) {
                    double tmp = options.mc3_temperatures[i];
                    if((tmp < 0.0) or (tmp > 1.0)) {
                        fprintf(stderr, "%s: temperatures must be 0.0 - 1.0 inclusive (%d%s value is %.3f)\n", 
                            argv[0], i+1, i==0 ? "st" : i==1 ? "nd" : i==2 ? "rd" : "th", tmp);
                        exit(EXIT_FAILURE);
                    }
                }

                break;
            
            case 'M':
                options.mc3 = true;
                break;
            */
            case 'a':
                options.affected_only = true;
                break;

            case 'e':
                options.elod = true;
                break;

            case 'w':
                if(not str2float(options.elod_marker_separation, optarg)) {
                    fprintf(stderr, "%s: option '-w' requires a float as an argument ('%s' given)\n", argv[0], optarg);
                    exit(EXIT_FAILURE);
                }
                if(options.elod_marker_separation <= 0.0) {
                    fprintf(stderr, "%s: ELOD marker separation must be greater than zero ('%f' given)\n", argv[0], options.elod_marker_separation);
                    exit(EXIT_FAILURE);
                }
                break;

            case 'k':
                if(not csv2vec(options.elod_penetrance, optarg)) {
                    fprintf(stderr, "%s: penetraces must be comma delimited floats from 0.0 - 1.0 inclusive, e.g.: 0.0,0.0,1.0 ('%s' given)\n", argv[0], optarg);
                    exit(EXIT_FAILURE);
                }
                
                if(options.elod_penetrance.size() != 3) {
                    fprintf(stderr, "%s: penetrance requires 3 floats, found %d!\n", argv[0], int(options.elod_penetrance.size()));
                    exit(EXIT_FAILURE);
                }

                for(unsigned i = 0; i < options.elod_penetrance.size(); ++i) {
                    double tmp = options.elod_penetrance[i];
                    if((tmp < 0.0) or (tmp > 1.0)) {
                        fprintf(stderr, "%s: temperatures must be 0.0 - 1.0 inclusive (%d%s value is %.3f)\n", 
                            argv[0], i+1, i==0 ? "st" : i==1 ? "nd" : i==2 ? "rd" : "th", tmp);
                        exit(EXIT_FAILURE);
                    }
                }

                break;

            case 'f':
                if(not str2float(options.elod_frequency, optarg)) {
                    fprintf(stderr, "%s: option '-f' requires a float as an argument ('%s' given)\n", argv[0], optarg);
                    exit(EXIT_FAILURE);
                }
                if(options.elod_frequency <= 0.0) {
                    fprintf(stderr, "%s: ELOD trait frequency must be greater than zero ('%f' given)\n", argv[0], options.elod_frequency);
                    exit(EXIT_FAILURE);
                }
                break;


            case 'u':
                if(not str2int(options.elod_replicates, optarg)) {
                    fprintf(stderr, "%s: option '-u' requires an int as an argument ('%s' given)\n", argv[0], optarg);
                    exit(EXIT_FAILURE);
                }
                if(options.elod_replicates < 1) {
                    fprintf(stderr, "%s: number of ELOD replicates must be greater than zero ('%d' given)\n", argv[0], options.elod_replicates);
                    exit(EXIT_FAILURE);
                }
                break;

            case 'X':
                options.sex_linked = true;
                break;

            case 'R':
                if(not str2int(options.mcmc_runs, optarg)) {
                    fprintf(stderr, "%s: option '-R' requires an int as an argument ('%s' given)\n", argv[0], optarg);
                    exit(EXIT_FAILURE);
                }
                if(options.mcmc_runs < 1) {
                    fprintf(stderr, "%s: number of MCMC runs must be greater than zero ('%d' given)\n", argv[0], options.mcmc_runs);
                    exit(EXIT_FAILURE);
                }
                break;

            case 'C':
                options.coda_logging = true;
                break;

            case 'P':
                options.coda_prefix = string(optarg);
                break;

            case ':':
                fprintf(stderr, "%s: option '-%c' requires an argument\n", 
                        argv[0], optopt);
                exit(EXIT_FAILURE);
                
			case '?':
			default:
                fprintf(stderr, "%s: option '-%c' invalid, ignoring...\n", 
						argv[0], optopt);
                break;
		}
	}

    //fprintf(stderr, "%s", options.debug_string().c_str());
	
    if(options.elod and pedfile != NULL)
        return;

    if(options.use_gpu and options.sex_linked) {
        fprintf(stderr, "Error: we do not current support sex-linked analysis on GPU\n");
        exit(EXIT_FAILURE);
    }

    if((int(options.mc3_temperatures.size()) > 0) and (int(options.mc3_temperatures.size()) != options.mc3_number_of_chains)) {
        fprintf(stderr, "Error: %d temperature%s %s set for %d chain%s\n", 
                int(options.mc3_temperatures.size()), 
                int(options.mc3_temperatures.size()) == 1 ? "" : "s", 
                int(options.mc3_temperatures.size()) == 1 ? "was" : "were",
                options.mc3_number_of_chains, 
                options.mc3_number_of_chains == 1 ? "" : "s");
        exit(EXIT_FAILURE);
    }

    if((int(options.mc3_temperatures.size()) > 0) and (options.mc3_temperatures[0] != 1.0)) {
        fprintf(stderr, "Warning: no sampling will happen if the temperature of the first chain is not 1.0!\n");
    }

	if((mapfile == NULL) or (pedfile == NULL) or (datfile == NULL)) {
	    fprintf(stderr, "Error: you must specify at least a pedigree file, a map file and a data file in LINKAGE format\n");
        exit(EXIT_FAILURE);
	}
}

void _set_runtime_parameters() {
    printf("setting %d thread%s\n", options.thread_count, options.thread_count == 1 ? "" : "s");
    set_num_threads(options.thread_count);
}

int linkage_analysis() {
    LinkageProgram lp(pedfile, mapfile, datfile, outfile, options);
    
    return lp.run() ? EXIT_SUCCESS : EXIT_FAILURE;
}

int elod_analysis() {
    Elod e(pedfile, options);

    double elod = e.run();

    fprintf(stderr, "\nELOD = %.3f\n", elod);
    
    return EXIT_SUCCESS;
}

int haplotype_analysis() {
/*
    HaplotypeProgram hp(pedfile, mapfile, datfile, verbose);

    return hp.run() ? EXIT_SUCCESS : EXIT_FAILURE;
*/

    fprintf(stderr, "Not currently supported...\n");

	return EXIT_FAILURE;
}

int testing_mode() {
    //TestProgram tp(pedfile, mapfile, datfile, outfile, mcmc_iterations, verbose);
    //return tp.run() ? EXIT_SUCCESS : EXIT_FAILURE;
    return EXIT_FAILURE;
}

int main(int argc, char **argv) {
    
#ifdef __linux__
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
#endif
    
    
	_handle_args(argc, argv);
	
	_set_runtime_parameters();
	
    if(options.elod)
        return elod_analysis();
	
	return linkage_analysis();
    
    /*
    switch(analysis) {
        case LINKAGE :
            return linkage_analysis();
            
        case HAPLOTYPE :
            return haplotype_analysis();
            
        case TESTING :
            return testing_mode();
        
        default :
            break;
    }

	return EXIT_FAILURE;
	*/
}

