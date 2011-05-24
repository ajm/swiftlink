using namespace std;

#include <cstdio>
#include <vector>
#include <cmath>

#include "pedigree.h"
#include "pedigree_parser.h"
#include "genetic_map.h"
#include "map_parser.h"
#include "disease_model.h"
#include "linkage_parser.h"
#include "peeler.h"


int main(int argc, char **argv) {
    
	if(argc != 4) {
		fprintf(stderr, "Usage: %s <ped> <map> <dat>\n", argv[0]);
		return EXIT_FAILURE;
	}

    DiseaseModel dm;
    
    GeneticMap map;
    MapParser mp(argv[2], map);
    
    if(!mp.parse()) {
        fprintf(stderr, "error parsing %s...\n", argv[2]);
		return EXIT_FAILURE;
    }
    
    LinkageParser lp(argv[3], map, dm);
	if(! lp.parse()) {
		fprintf(stderr, "errors parsing linkage file, exiting...\n");
		return EXIT_FAILURE;
	}
    
    if(! map.sanity_check()) {
		fprintf(stderr, "map data failed sanity check\n");
		return EXIT_FAILURE;
	}
    
    vector<Pedigree> v;
	PedigreeParser p(argv[1], v, dm);
	
	if(!p.parse()) {
		fprintf(stderr, "error parsing %s...\n", argv[1]);
		return EXIT_FAILURE;
	}
    

    Peeler test(&v[0], &map);
    double trait_prob = test.get_trait_prob();
    
    printf("P(T) = %e\n", trait_prob / log(10));
    
    return EXIT_SUCCESS;
}
