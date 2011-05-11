using namespace std;

#include <cstdio>
#include <vector>
#include <cmath>

#include "peel_sequence_generator.h"
#include "peeling.h"
#include "peeler.h"
#include "genetic_map.h"
#include "pedigree.h"
#include "descent_graph.h"
#include "rfunction.h"


Peeler::Peeler(Pedigree* p, GeneticMap* g) 
    : ped(p), map(g), lod(p, g) {
    
    PeelSequenceGenerator psg(*ped);
    psg.build_peel_order();
    psg.print();

    vector<PeelOperation>& ops = psg.get_peel_order();
    
    for(vector<PeelOperation>::size_type i = 0; i < ops.size(); ++i) {
        Rfunction rf(ops[i], ped, map, 4);
        rfunctions.push_back(rf);
    }
    
    lod.set_trait_prob(calc_trait_prob());
}

double Peeler::calc_trait_prob() {
    return this->peel(NULL, 0);
}

bool Peeler::process(DescentGraph* dg) {
    // minus 1 because we want to look between markers
    // m-t-m-t-m-t-m where m is a marker and t is a trait location
    for(unsigned i = 0; i < map->num_markers() - 1; ++i) {
        lod.add(i, this->peel(dg, i) - dg->trans_prob()); // this is all in base e
    }
    
    return true;
}

double Peeler::peel(DescentGraph* dg, unsigned locus) {
    PeelMatrix* last = NULL;

    for(unsigned i = 0; i < rfunctions.size(); ++i) {
        Rfunction& rf = rfunctions[i];
        
        if(not rf.evaluate(last, dg, locus)) {
            return false;
        }
        
        last = rf.get_matrix();
    }
    
    return log(last->get_result()); // TODO : convert all peel code to log likelihood?
}

double Peeler::get(unsigned locus) {
    return lod.get(locus);
}

