using namespace std;

#include <cstdio>
#include <cstdlib>
//#include <ctime>
#include <vector>

#include "linkage_program.h"

#include "genetic_map.h"
#include "pedigree.h"
#include "peel_sequence_generator.h"
#include "peeling.h"
#include "peel_matrix.h"
#include "rfunction.h"
//#include "descentgraph.h"


bool LinkageProgram::run() {
	//SA_DescentGraph* ans;

	//srand(time(NULL));

//	for(unsigned int i = 0; i < pedigrees.size(); ++i) {
//        fprintf(stderr, "processing pedigree %d (id = %s)\n", i, pedigrees[i].get_id().c_str());
        
		//SimulatedAnnealing sa(pedigrees[i], &map);
		//ans = sa.optimise(1600 * pedigrees[i]->num_members() * pedigrees[i]->num_markers() * 20 * 2);

		// ???
		
		//delete ans;
//	}

    // for each pedigree
	for(vector<Pedigree>::size_type i = 0; i < pedigrees.size(); ++i) {

        fprintf(stderr, "processing pedigree %d (id = %s)\n", i, pedigrees[i].get_id().c_str());

        PeelSequenceGenerator psg(pedigrees[i]);
        psg.build_peel_order();
        psg.print();
        
        // setup r-functions
        vector<PeelOperation>& ops = psg.get_peel_order();
        
        
        vector<Rfunction> rfunctions;
        for(vector<PeelOperation>::size_type j = 0; j < ops.size(); ++j) {
            Rfunction rf(ops[j], &pedigrees[i], 4);
            rfunctions.push_back(rf);
        }

        // perform the peel for every locus
        PeelMatrix* last = NULL;
        for(vector<Rfunction>::size_type j = 0; j < rfunctions.size(); ++j) {
            Rfunction& rf = rfunctions[j];

            fprintf(stderr, "rfunction %d\n", int(j));
            
            if(not rf.evaluate(last)) {
                fprintf(stderr, "bad R function\n");
                return 1;
            }
            
            last = rf.get_matrix();
        }
	}

	return true;
}

