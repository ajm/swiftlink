using namespace std;

#include <cstdio>
#include <vector>

#include "parser.h"
#include "pedigree.h"
#include "pedfileparser.h"
#include "peel_sequence_generator.h"
#include "diseasemodel.h"
//#include "rfunction.h"


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

    // for each pedigree
	for(vector<Pedigree>::size_type i = 0; i < v.size(); ++i) {
        PeelSequenceGenerator psg(v[i]);
        psg.build_peel_order();
        
        // setup r-functions
        vector<PeelOperation>& ops = psg.get_peel_order();

        // NONE OF THIS EVEN COMPILES  ;-P
/*        
        vector<Rfunction> rfunctions;
        for(vector<PeelOperation>::size_type j = 0; j < ops.size(); ++j) {
            Rfunction rf(ops[j], &v[i], 4);
            rfunctions.push_back(rf);
        }

        // perform the peel
        PeelMatrix* last = NULL;
        for(vector<Rfunction>::size_type j = 0; j < rfunctions.size(); ++j) {
            Rfunction& rf = rfunctions[j];
            
            rf.evaluate(last);
            
            last = rf.get_peel_matrix();
        }
*/
	}

	return 0;
}

