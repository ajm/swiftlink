using namespace std;

#include <cstdio>
#include <vector>
#include <cmath>

#include "peel_sequence_generator.h"
#include "peeling.h"
#include "peeler.h"
#include "genetic_map.h"
#include "pedigree.h"
#include "simwalk_descent_graph.h"
#include "rfunction.h"


Peeler::Peeler(Pedigree* p, GeneticMap* g) : ped(p), map(g), lod(p, g) {
    
    PeelSequenceGenerator psg(*ped);
    psg.build_peel_order();
    psg.print();

    vector<PeelOperation>& ops = psg.get_peel_order();
    
    for(vector<PeelOperation>::size_type i = 0; i < ops.size(); ++i) {
        Rfunction rf(ops[i], ped, map, 4);
        rfunctions.push_back(rf);
    }
}

// XXX put all peeling sequence calculation + evaluation code 
// in own module, this is not the final form as it need to peel
// over many descent graphs + all loci not just this one, but this is just
// placeholder to get started...
bool Peeler::peel(SimwalkDescentGraph* sdg) {
        
    printf("\n");
    sdg->print();

    // minus 1 because we want to look between markers
    // m-t-m-t-m-t-m where m is a marker and t is a trait location
    for(unsigned i = 0; i < map->num_markers() - 1; ++i) {
        PeelMatrix* last = NULL;

        for(unsigned j = 0; j < rfunctions.size(); ++j) {
            Rfunction& rf = rfunctions[j];
            
            if(not rf.evaluate(last, sdg, i)) {
                return false;
            }
                        
            last = rf.get_matrix();
            
            printf("\nRfunction %d\n\n", int(j));
            last->print();
        }
        
        lod.add(i, last->get_result(), sdg->trans_prob() / log(10));
        
        printf("locus = %d, prob = %f (P(G^) = %f, Trans(G^) = %f)\n", 
               i, 
               log(last->get_result()), 
               sdg->get_prob(), 
               sdg->trans_prob()
            );

        //ajm
        exit(-1);
    }
    
    return true;
}

void Peeler::print() {
    lod.print();
}
