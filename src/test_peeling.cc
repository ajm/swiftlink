using namespace std;

#include <cstdio>
#include <vector>

#include "parser.h"
#include "pedigree.h"
#include "pedfileparser.h"
#include "peel_sequence_generator.h"
#include "diseasemodel.h"


// XXX
// XXX this is all very clumsy
// XXX
int main(int argc, char **argv) {

	if(argc != 2) {
		fprintf(stderr, "Usage: %s <pedfile>\n", argv[0]);
		return 1;
	}
    
    DiseaseModel dm;
    dm.set(SIMPLE_AUTOSOMAL_RECESSIVE);

	vector<Pedigree> v;
	PedigreeParser p(argv[1], v, dm);
	
	if(!p.parse()) {
		fprintf(stderr, "parsing error\n");
		return 1;
	}
/*
    // for each pedigree in pedfile
	for(unsigned int i = 0; i < v.size(); ++i) {
        // find a suitable peel order
        PeelSequenceGenerator psg(v[i]);
        psg.build_peel_order();
        psg.print();
        
        // perform the peel
        // ???
	}
*/
	return 0;
}

