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
    
    bool sane = true;
	vector<Pedigree*> v;
	PedigreeParser p(argv[1], &v);
	
	if(!p.parse()) {
		fprintf(stderr, "parsing error\n");
		return 1;
	}

    DiseaseModel dm;
    dm.set_freq(1 / 1000.0);
    dm.set_sexlinked(false);
	dm.set_penetrance(1.0, 0);
    dm.set_penetrance(1.0, 1);
    dm.set_penetrance(1.0, 2);
	
    // for each pedigree in pedfile
	for(unsigned int i = 0; i < v.size(); ++i) {
        // print if sane
		if(v[i]->sanity_check())
			v[i]->print();
		else {
			sane = false;
			fprintf(stderr, "sanity error\n");
		}

        // set person local probabilities
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
	
	return sane ? 0 : 1;
}

