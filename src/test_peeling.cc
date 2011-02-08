using namespace std;

#include <cstdio>
#include <vector>

#include "parser.h"
#include "pedigree.h"
#include "pedfileparser.h"
#include "peel_sequence_generator.h"

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
	
	for(unsigned int i = 0; i < v.size(); ++i) {
		if(v[i]->sanity_check()) {
			v[i]->print();
		}
		else {
			sane = false;
			fprintf(stderr, "sanity error\n");
		}

        // perform peel
        PeelSequenceGenerator psg(v[i]);
        psg.build_peel_order();
        psg.print();
        

		delete v[i]; // just to shut up valgrind
	}
	
	return sane ? 0 : 1;
}

