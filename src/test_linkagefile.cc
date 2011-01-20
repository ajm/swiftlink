using namespace std;

#include <cstdio>

#include "geneticmap.h"
#include "diseasemodel.h"
#include "mapfileparser.h"
#include "linkagefileparser.h"

int main(int argc, char **argv) {
	
    if(argc != 3) {
		fprintf(stderr, "Usage: %s <mapfile> <linkagefile>\n", argv[0]);
		return 1;
	}
	
    GeneticMap m;
    MapParser p(argv[1], &m);
	
    if(!p.parse()) {
        fprintf(stderr, "map file parsing error\n");
        return 1;
    }
	
    printf("read %d markers\n", int(m.num_markers()));
	
	DiseaseModel d;
	LinkageParser lp(argv[2], &m, &d);
	
	if(!lp.parse()) {
		fprintf(stderr, "linkage file parsing error\n");
        return 1;
	}
	
	d.print();
	m.print();
	
    return 0;
}

