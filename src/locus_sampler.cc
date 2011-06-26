using namespace std;

#include <vector>

#include "descent_graph_diff.h"
#include "descent_graph.h"
#include "peel_sequence_generator.h"
#include "peeling.h"
#include "pedigree.h"
#include "genetic_map.h"
#include "rfunction.h"
#include "rfunction_builder.h"
#include "sampler_rfunction.h"
#include "locus_sampler.h"


LocusSampler::LocusSampler(Pedigree* ped, GeneticMap* map, DescentGraph* dg) :
    Sampler(ped, map, dg), 
    rfunctions() {
    
    PeelSequenceGenerator psg(ped);
    psg.build_peel_order();

    vector<PeelOperation>& ops = psg.get_peel_order();
    
    RfunctionBuilder<SamplerRfunction> build(ped, map, rfunctions);
        
    for(unsigned i = 0; i < ops.size(); ++i) {
        rfunctions.push_back(build.createRfunction(ops[i]));
    }
}

void LocusSampler::step(DescentGraphDiff& dgd) {
    unsigned locus = get_random_locus();
    
    // forward peel
    for(unsigned i = 0; i < rfunctions.size(); ++i) {
        SamplerRfunction* rf = rfunctions[i];
        rf->evaluate(dg, locus, 0.0);
    }
/*    
    // reverse peel
    Map<unsigned, enum phased_trait> assignment;
    for(unsigned i = rfunctions.size(); i >= 0; --i) {
        SamplerRfunction* rf = rfunctions[i];
        rf->sample(assignment);
    }
*/
    // return a DescentGraphDiff object
    
}

