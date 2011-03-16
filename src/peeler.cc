using namespace std;

#include <cstdio>
#include <vector>

#include "peel_sequence_generator.h"
#include "peeling.h"
#include "peeler.h"
#include "genetic_map.h"
#include "pedigree.h"
#include "simwalk_descent_graph.h"
#include "progress.h"
#include "rfunction.h"


Peeler::Peeler(Pedigree* p, GeneticMap* g) : ped(p), map(g) {
    
    PeelSequenceGenerator psg(*ped);
    psg.build_peel_order();
    psg.print();

    //exit(0);

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
    
    Progress prog("Peeling:", rfunctions.size() * map->num_markers());

    // minus 1 because we want to look between markers
    // m-t-m-t-m-t-m where m is a marker and t is a trait location
    for(unsigned int i = 0; i < map->num_markers() - 1; ++i) {
        
        fprintf(stderr, "\n\nPeeling locus %d\n\n", int(i));
        
        PeelMatrix* last = NULL;

        for(vector<Rfunction>::size_type j = 0; j < rfunctions.size(); ++j) {
            Rfunction& rf = rfunctions[j];

            // XXX rf is getting reused over and over again
            // 'clear'? contains or set a dirty bit somewhere?
            
            if(not rf.evaluate(last, sdg, i)) {
                prog.error("bad R function");
                return false;
            }
            
            prog.increment();
            
            last = rf.get_matrix();
        }

        // XXX print the likelihood?
    }

    prog.finish();

    return true;
}

