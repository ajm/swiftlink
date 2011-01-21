using namespace std;

#include <cstdio>
#include <vector>

#include "parser.h"
#include "pedfileparser.h"
#include "pedigree.h"
#include "peel.h"

int main(int argc, char **argv) {

	if(argc != 2) {
		fprintf(stderr, "Usage: %s <pedfile>\n", argv[0]);
		return 1;
	}

    bool sane = true;
	vector<Pedigree*> pv;
	PedigreeParser p(argv[1], &pv);
	
	if(!p.parse()) {
		fprintf(stderr, "parsing error\n");
		return 1;
	}
	
	//printf("read in %u families\n", pv.size());
	
	for(unsigned int i = 0; i < pv.size(); ++i) {
		if(pv[i]->sanity_check()) {
			pv[i]->print();
		}
		else {
			sane = false;
			fprintf(stderr, "sanity error\n");
		}

        // perform peel

		delete pv[i]; // just to shut up valgrind
	}
	
	return sane ? 0 : 1;
}

