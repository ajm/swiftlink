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
#include "rfunction_builder.h"


Peeler::Peeler(Pedigree* p, GeneticMap* g) : 
    ped(p), 
    map(g), 
    rfunctions(),
    lod(p, g),
    trait_prob(0.0) {
    
    PeelSequenceGenerator psg(ped);
    psg.build_peel_order();
    psg.print();

    vector<PeelOperation>& ops = psg.get_peel_order();
        
    RfunctionBuilder<TraitRfunction> build(ped, map, rfunctions);
    
    for(unsigned i = 0; i < ops.size(); ++i) {
        rfunctions.push_back(build.createRfunction(ops[i]));
    }
    
    trait_prob = calc_trait_prob();
    
    printf("P(T) = %e\n", trait_prob / log(10));
    printf("PEEL SEQUENCE SCORE (lower is better) : %d\n", int(psg.score_peel_sequence()));
    
    lod.set_trait_prob(trait_prob);
}

Peeler::Peeler(const Peeler& rhs) :
    ped(rhs.ped),
    map(rhs.map),
    rfunctions(),
    lod(rhs.lod),
    trait_prob(rhs.trait_prob) {

    copy_rfunctions(rhs);
}

Peeler& Peeler::operator=(const Peeler& rhs) {
    
    if(&rhs != this) {
        ped = rhs.ped;
        map = rhs.map;
        copy_rfunctions(rhs);
        lod = rhs.lod;
        trait_prob = rhs.trait_prob;
    }
    
    return *this;
}

Peeler::~Peeler() {
    for(unsigned i = 0; i < rfunctions.size(); ++i)
        delete rfunctions[i];
}

void Peeler::copy_rfunctions(const Peeler& rhs) {
    rfunctions.clear();
    for(unsigned i = 0; i < rhs.rfunctions.size(); ++i) {
        TraitRfunction* rf = new TraitRfunction(*(rhs.rfunctions[i]));
        rfunctions.push_back(rf);
    }
}

double Peeler::get_trait_prob() {
    return trait_prob;
}

double Peeler::calc_trait_prob() {
    return peel(NULL, 0);
}

bool Peeler::process(DescentGraph& dg) {
    // minus 1 because we want to look between markers
    // m-t-m-t-m-t-m where m is a marker and t is a trait location
    for(unsigned i = 0; i < map->num_markers() - 1; ++i) {
        lod.add(i, peel(&dg, i) - dg.trans_prob()); // this is all in base e
    }
    
    return true;
}

double Peeler::peel(DescentGraph* dg, unsigned locus) {

    for(unsigned i = 0; i < rfunctions.size(); ++i) {
        TraitRfunction* rf = rfunctions[i];
        rf->evaluate(dg, locus, 0.5); // TODO XXX <---- 0.5 is not actually used yet...
    }
    
    TraitRfunction* rf = rfunctions.back();
    
    return log(rf->get_result()); // TODO : convert all peel code to log likelihood?
}

// XXX change this so the peeler has the responsibility to
// print everything to a file
double Peeler::get(unsigned locus) {
    return lod.get(locus);
}

