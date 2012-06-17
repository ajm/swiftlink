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

#include <omp.h>

#include "types.h"
#include "defaults.h"
#include "linkage_program.h"
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
"\n"
"Runtime options:\n"
"  -c NUM,     --cores=NUM                 (default = %d)\n"
"  -g,         --gpu\n"
"\n"
"Misc:\n"
"  -q seqfile, --peelsequence=seqfile\n"
"  -v,         --verbose\n"
"  -h,         --help\n"
"\n", 
prog, 
DEFAULT_RESULTS_FILENAME, 
DEFAULT_MCMC_ITERATIONS,
DEFAULT_BURNIN_ITERATIONS,
DEFAULT_SEQUENTIALIMPUTATION_RUNS,
DEFAULT_SCORING_PERIOD,
DEFAULT_LSAMPLER_PROB,
DEFAULT_LODSCORES,
DEFAULT_THREAD_COUNT
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

void _handle_args(int argc, char **argv) {
	extern char *optarg;
    extern int optopt;
	int ch;
	
	static struct option long_options[] = 
	    {
	        {"verbose",             no_argument,        0,      'v'},
	        {"help",                no_argument,        0,      'h'},
	        {"iterations",          required_argument,  0,      'i'},
	        {"burnin",              required_argument,  0,      'b'},
	        {"sequentialimputation",required_argument,  0,      's'},
	        {"lsamplerprobability", required_argument,  0,      'l'},
	        {"scoringperiod",       required_argument,  0,      'x'},
	        {"gpu",                 no_argument,        0,      'g'},
	        {"cores",               required_argument,  0,      'c'},
	        {"pedigree",            required_argument,  0,      'p'},
	        {"map",                 required_argument,  0,      'm'},
	        {"dat",                 required_argument,  0,      'd'},
	        {"output",              required_argument,  0,      'o'},
	        {"peelsequence",        required_argument,  0,      'q'},
	        {"lodscores",           required_argument,  0,      'n'},
	        {0, 0, 0, 0}
	    };
    
    int option_index = 0;
    
	while ((ch = getopt_long(argc, argv, ":p:d:m:o:i:b:s:l:c:x:q:n:vhcg", long_options, &option_index)) != -1) {
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
                if(options.si_iterations < 0) {
                    fprintf(stderr, "%s: number of sequential imputation runs must be positive (%d given)\n", argv[0], options.si_iterations);
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
                options.use_gpu = true;
                break;
                
            case 'q':
                options.peelseq_filename = string(optarg);
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
	
	if((mapfile == NULL) or (pedfile == NULL) or (datfile == NULL)) {
	    fprintf(stderr, 
	    "\n"
	    "-----------------------------------------------------------\n"
	    "|                                                         |\n"
	    "|  error: at a minimum you must specify a pedigree file,  |\n"
	    "|         a map file and a linkage format dat file        |\n"
	    "|                                                         |\n"
	    "-----------------------------------------------------------\n"
	    "\n");
        _usage(argv[0]);
        exit(EXIT_FAILURE);
	}
}

void _set_runtime_parameters() {
    printf("setting %d threads\n", options.thread_count);
    omp_set_num_threads(options.thread_count);
}

int linkage_analysis() {
    LinkageProgram lp(pedfile, mapfile, datfile, outfile, options);
    
    return lp.run() ? EXIT_SUCCESS : EXIT_FAILURE;
}

int haplotype_analysis() {
/*
    HaplotypeProgram hp(pedfile, mapfile, datfile, verbose);

    return hp.run() ? EXIT_SUCCESS : EXIT_FAILURE;
*/

    fprintf(stderr, "Sorry chaps, not currently supported...\n");

	return EXIT_FAILURE;
}

int testing_mode() {
    //TestProgram tp(pedfile, mapfile, datfile, outfile, mcmc_iterations, verbose);
    //return tp.run() ? EXIT_SUCCESS : EXIT_FAILURE;
    return EXIT_FAILURE;
}

int main(int argc, char **argv) {
    
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
    
    
	_handle_args(argc, argv);
	
	_set_runtime_parameters();
	
	
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

