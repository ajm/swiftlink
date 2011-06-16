/*
 * the only thing this file does is handle command line arguments
 * and then call another main function
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#include "misc.h"
#include "linkage_program.h"
#include "haplotype_program.h"


enum analysistype {
    LINKAGE,
    HAPLOTYPE
};

char*   mapfile;
char*   pedfile;
char*   datfile;
bool    verbose;
enum analysistype analysis;


void _usage(char *prog) {
	fprintf(stderr,
			"Usage: %s [-hv] -p pedigreefile -m mapfile -d datfile\n"
			"\t-p <pedigree file>\n"
			"\t-m <map file>\n"
			"\t-d <data file>\n"
			"\t-a <linkage|haplotype> (default: linkage)\n"
			"\t-v verbose\n"
			"\t-h print usage information\n"
			"\n", 
			prog);
}

void _meow(void) {
	fprintf(stderr,
			" __________\n"
			"|          |\n"
			"|   meow   |\n"
			"|____   ___|\n"
			"      \\| _\n"
			"         \\`.|\\\n"
			"         /  ' `\n"
			"         )/' _/\n"
			"         `-'\"\n"
			);
	
	exit(EXIT_SUCCESS);
}

void _handle_args(int argc, char **argv) {
	extern char *optarg;
    extern int optopt;
	int ch;
	int bad = 0;
	
	mapfile = pedfile = datfile = NULL;
	verbose = false;
	analysis = LINKAGE;
	
	while ((ch = getopt(argc, argv, ":p:d:m:a:vhc")) != -1) {
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
				verbose = true;
				break;
				
		    case 'a':
		        if(strcmp(optarg, "linkage") == 0) {
		            analysis = LINKAGE;
		        }
		        else if(strcmp(optarg, "haplotype") == 0) {
		            analysis = HAPLOTYPE;
		        }
		        else {
		            fprintf(stderr, "%s: option '-a' can only accept 'linkage' or "
		                "'haplotype' as arguments ('%s' given)\n",
		                argv[0], optarg);
		            exit(-1);
		        }
		        break;
		        
			case 'h':
				_usage(argv[0]);
				exit(EXIT_SUCCESS);
				
			case 'c':
				_meow();
				
            case ':':
                fprintf(stderr, "%s: option '-%c' requires an argument\n", 
                        argv[0], optopt);
                break;
                
			case '?':
			default:
                fprintf(stderr, "%s: option '-%c' invalid, ignoring...\n", 
						argv[0], optopt);
                break;
		}
	}
	
	if(mapfile == NULL) {
		fprintf(stderr, "error: you must specify a map file\n");
		bad = 1;
	}
	
	if(pedfile == NULL) {
		fprintf(stderr, "error: you must specify a pedigree file\n");
		bad = 1;
	}
	
	if(datfile == NULL) {
		fprintf(stderr, "error: you must specify a data file\n");
		bad = 1;
	}
	
	if(bad)
		exit(EXIT_FAILURE);
}

int linkage_analysis() {
    LinkageProgram lp(pedfile, mapfile, datfile, verbose);
    
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

int main(int argc, char **argv) {
	_handle_args(argc, argv);
    
    switch(analysis) {
        case LINKAGE :
            return linkage_analysis();
            
        case HAPLOTYPE :
            return haplotype_analysis();
        
        default :
            break;
    }

	return EXIT_FAILURE;
}

