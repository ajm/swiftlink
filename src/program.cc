using namespace std;

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

    if(! map.sanity_check()) {
		fprintf(stderr, "error: map data failed sanity check...\n");
		return false;
	}

    // read in the pedigree file
	PedigreeParser pp(pedfile, pedigrees, dm, map);
	if(! pp.parse()) {
		//fprintf(stderr, "errors parsing pedigree file, exiting...\n");
		return false;
	}
    
	return true;
}

