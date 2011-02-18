#include <cstdio>

#include "genetic_map.h"


bool GeneticMap::sanity_check() {
    bool sane = (map.size() == (thetas.size() - 1)) and \
           (map.size() == (inverse_thetas.size() - 1));

    if(not sane) {
        fprintf(stderr, 
            "error in map data: number of markers = %d, number of thetas = %d\n", 
            map.size(), thetas.size());
    }
    
    return sane;
}

void GeneticMap::print() {
    printf("GeneticMap: %d loci\n", map.size());
	
	for(unsigned int i = 0; i < map.size(); ++i)
		map[i].print();
		
	printf("\n");
	printf("  thetas:\n");
	for(unsigned int i = 0; i < thetas.size(); ++i)
		printf("\t%f\n", exp(thetas[i]));

	printf("\n");
}

