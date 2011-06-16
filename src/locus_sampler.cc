using namespace std;

#include <vector>

#include "descent_graph_diff.h"
#include "descent_graph.h"
#include "peel_sequence_generator.h"
#include "rfunction.h"
#include "peeling.h"
#include "pedigree.h"
#include "genetic_map.h"


LocusSampler::LocusSampler(Pedigree* ped, DescentGraph* dg) :
    Sampler(ped, dg) {
    
    PeelSequenceGenerator psg(ped);
    psg.build_peel_order();

    vector<PeelOperation>& ops = psg.get_peel_order();
    
    for(vector<PeelOperation>::size_type i = 0; i < ops.size(); ++i) {
        Rfunction* rf = new Rfunction(ops[i], ped, map, 4, rfunctions, i);
        rfunctions.push_back(rf);
    }
}

void LocusSampler::step(DescentGraphDiff& dgd) {
    unsigned locus = get_random_locus();
    
    // forward peel
    for(unsigned i = 0; i < rfunctions.size(); ++i) {
        Rfunction* rf = rfunctions[i];
        rf->evaluate(dg, locus, 0.0);
    }
    
    // reverse peel
    
    // return a DescentGraphDiff object
    
}

