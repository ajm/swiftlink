#include <cstdio>
#include <vector>

#include "types.h"
#include "program.h"

#include "linkage_parser.h"
#include "pedigree_parser.h"
#include "map_parser.h"

#include "disease_model.h"
#include "genetic_map.h"
#include "pedigree.h"

using namespace std;


bool Program::read_and_check_input() {

    // read in the map file
    MapParser mp(mapfile, map);
	if(! mp.parse()) {
		//fprintf(stderr, "errors parsing map file, exiting...\n");
		return false;
	}
	
    // read 'linkage' file format config file
    LinkageParser lp(datfile, map, dm);
	if(! lp.parse()) {
		//fprintf(stderr, "errors parsing linkage file, exiting...\n");
		return false;
	}

    if(options.verbose) {
        fprintf(stderr, "%s\n", dm.debug_string().c_str());
    }

    if(options.sex_linked) {
        dm.set_sexlinked(true);
    }

    if(! map.sanity_check()) {
		fprintf(stderr, "Error: map data failed sanity check...\n");
		return false;
	}

    // read in the pedigree file
	PedigreeParser pp(pedfile, pedigrees, dm, map);
	if(! pp.parse()) {
		//fprintf(stderr, "errors parsing pedigree file, exiting...\n");
		return false;
	}

    for(unsigned int i = 0; i < pedigrees.size(); ++i) {
        if(map.num_markers() != pedigrees[i].num_markers()) {
            fprintf(stderr, 
                    "Error: different number of markers in \"%s\" pedigree (%d) versus map (%d)\nExiting...\n", 
                    pedigrees[i].get_id().c_str(), pedigrees[i].num_markers(), map.num_markers());
            exit(EXIT_FAILURE);
        }
    }

	return true;
}

