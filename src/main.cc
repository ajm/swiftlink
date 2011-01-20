/*
 * the only thing this file does is handle command line arguments
 * and then call another main function
 */

#include <cstdio>
#include <cstdlib>
#include <unistd.h>

#include "linkageprogram.h"


char* mapfile;
char* pedfile;
char* datfile;
bool verbose;

void _usage(char *prog) {
	fprintf(stderr,
			"usage: %s [-hv] -p pedigreefile -m mapfile -d datfile\n"
			"\t-p <pedigree file>\n"
			"\t-m <map file>\n"
			"\t-d <data file>\n"
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

static void _handle_args(int argc, char **argv) {
	extern char *optarg;
    extern int optopt;
	int ch;
	int bad = 0;
	
	mapfile = pedfile = datfile = NULL;
	verbose = false;
	
	while ((ch = getopt(argc, argv, ":p:d:m:vhc")) != -1) {
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
	
	if(bad) {
		//_usage(argv[0]);
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char **argv) {
	_handle_args(argc, argv);
    
    LinkageProgram lp(pedfile, mapfile, datfile);
	if(lp.run()) {
		return EXIT_SUCCESS;
	}

	return EXIT_FAILURE;
}

