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
#include "trait_rfunction.h"


Peeler::Peeler(Pedigree* p, GeneticMap* g, PeelSequenceGenerator& psg) : 
    ped(p), 
    map(g), 
    rfunctions(),
    lod(p, g),
    trait_prob(0.0) {
    
    vector<PeelOperation>& ops = psg.get_peel_order();
    
    rfunctions.reserve(ops.size()); // need to do this otherwise pointers may not work later...
    
    for(unsigned i = 0; i < ops.size(); ++i) {
        Rfunction* prev1 = ops[i].get_previous_op1() == -1 ? NULL : &rfunctions[ops[i].get_previous_op1()];
        Rfunction* prev2 = ops[i].get_previous_op2() == -1 ? NULL : &rfunctions[ops[i].get_previous_op2()];
        
        rfunctions.push_back(TraitRfunction(ped, map, 0, &(ops[i]), prev1, prev2));
    }
    
    trait_prob = calc_trait_prob();
    
    printf("P(T) = %e\n", trait_prob / log(10));
    printf("PEEL SEQUENCE SCORE (lower is better) : %d\n", int(psg.score_peel_sequence()));
    
    lod.set_trait_prob(trait_prob);
}

Peeler::Peeler(const Peeler& rhs) :
    ped(rhs.ped),
    map(rhs.map),
    rfunctions(rhs.rfunctions),
    lod(rhs.lod),
    trait_prob(rhs.trait_prob) {

    //copy_rfunctions(rhs);
}

Peeler& Peeler::operator=(const Peeler& rhs) {
    
    if(&rhs != this) {
        ped = rhs.ped;
        map = rhs.map;
        rfunctions = rhs.rfunctions;
        lod = rhs.lod;
        trait_prob = rhs.trait_prob;

        //kill_rfunctions();
        //copy_rfunctions(rhs);
    }
    
    return *this;
}

Peeler::~Peeler() {}

double Peeler::get_trait_prob() {
    return trait_prob;
}

double Peeler::calc_trait_prob() {
    return peel(NULL, 0);
}

void Peeler::process(DescentGraph& dg) {
    // minus 1 because we want to look between markers
    // m-t-m-t-m-t-m where m is a marker and t is a trait location
    for(unsigned i = 0; i < map->num_markers() - 1; ++i) {
    
        //printf("%.4f\n", peel(&dg, i));
    
        lod.add(i, peel(&dg, i) - dg.get_recombination_prob(i, false) - dg.get_marker_transmission());
    }
}

double Peeler::peel(DescentGraph* dg, unsigned locus) {

    for(unsigned i = 0; i < rfunctions.size(); ++i) {
        //TraitRfunction* rf = rfunctions[i];
        //rf->evaluate(dg, 0.5); // TODO XXX <---- 0.5 is not actually used yet, always assumed to be 0.5
        
        rfunctions[i].set_locus(locus);
        rfunctions[i].evaluate(dg, 0.5);
    }
    
    TraitRfunction& rf = rfunctions.back();
    
    return log(rf.get_result());
}

double Peeler::get(unsigned locus) {
    return lod.get(locus);
}

