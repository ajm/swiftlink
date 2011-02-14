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
    
	vector<Pedigree> v;
	PedigreeParser p(argv[1], v);
	
	if(!p.parse()) {
		fprintf(stderr, "parsing error\n");
		return 1;
	}

/*
    DiseaseModel dm;
    dm.set(SIMPLE_AUTOSOMAL_RECESSIVE);
	
    // for each pedigree in pedfile
	for(unsigned int i = 0; i < v.size(); ++i) {
        // set person local probabilities
        // XXX this should be done in the pedigree object
        // the disease model should be passed to the PedigreeParser object
        for(unsigned int j = 0; j < v[i]->num_members(); ++j) {
            Person* p = v[i]->get_by_index(j);
            p->init_probs(dm);
        }

        // find a suitable peel order
        PeelSequenceGenerator psg(v[i]);
        psg.build_peel_order();
        psg.print();

        // perform the peel
        

		delete v[i]; // just to shut up valgrind
	}
*/	
	return 0;
}

