using namespace std;

#include <cstdio>
#include <cstdlib>
//#include <ctime>
#include <vector>

#include "linkage_program.h"

#include "geneticmap.h"
#include "pedigree.h"
//#include "descentgraph.h"


bool LinkageProgram::run() {
	//SA_DescentGraph* ans;

	//srand(time(NULL));

	for(unsigned int i = 0; i < pedigrees.size(); ++i) {
        fprintf(stderr, "processing pedigree %d (id = %s)\n", i, pedigrees[i].get_id().c_str());
        
		//SimulatedAnnealing sa(pedigrees[i], &map);
		//ans = sa.optimise(1600 * pedigrees[i]->num_members() * pedigrees[i]->num_markers() * 20 * 2);

		// ???
		
		//delete ans;
	}

	return true;
}

