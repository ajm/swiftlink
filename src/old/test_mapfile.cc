using namespace std;

#include <cstdio>

#include "geneticmap.h"
#include "mapfileparser.h"

int main(int argc, char **argv) {
    
    if(argc != 2) {
		fprintf(stderr, "Usage: %s <mapfile>\n", argv[0]);
		return 1;
	}
    
    GeneticMap m;
    MapParser p(argv[1], &m);
    
    if(!p.parse()) {
        fprintf(stderr, "parsing error\n");
        return 1;
    }
    
    printf("read %d markers\n", int(m.num_markers()));
    
    return 0;
}

