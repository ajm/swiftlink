using namespace std;

#include <cstdio>
#include <vector>

#include "parser.h"
#include "pedfileparser.h"
#include "pedigree.h"
#include "peeling.h"

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
	
	for(unsigned int i = 0; i < pv.size(); ++i) {
		if(v[i]->sanity_check()) {
			v[i]->print();
		}
		else {
			sane = false;
			fprintf(stderr, "sanity error\n");
		}

        // perform peel
        PedigreePeeler pp(v[i]);
        pp.build_peel_order();
        pp.print();
        

		delete v[i]; // just to shut up valgrind
	}
	
	return sane ? 0 : 1;
}

