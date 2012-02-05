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


Peeler::Peeler(Pedigree* p, GeneticMap* g, PeelSequenceGenerator* psg, unsigned int locus) : 
    ped(p), 
    map(g), 
    rfunctions(),
    trait_prob(0.0),
    lod_score(0.0),
    initialised(false),
    locus(locus),
    count(0) {
    
    vector<PeelOperation>& ops = psg->get_peel_order();
    
    rfunctions.reserve(ops.size()); // need to do this otherwise pointers may not work later...
    
    for(unsigned i = 0; i < ops.size(); ++i) {
        Rfunction* prev1 = ops[i].get_previous_op1() == -1 ? NULL : &rfunctions[ops[i].get_previous_op1()];
        Rfunction* prev2 = ops[i].get_previous_op2() == -1 ? NULL : &rfunctions[ops[i].get_previous_op2()];
        
        rfunctions.push_back(TraitRfunction(ped, map, 0, &(ops[i]), prev1, prev2));
    }
    
    trait_prob = calc_trait_prob();
    
    for(unsigned i = 0; i < rfunctions.size(); ++i) {
        rfunctions[i].set_locus(locus);
    }
    
    //printf("P(T) = %e\n", trait_prob / log(10));
    //printf("PEEL SEQUENCE SCORE (lower is better) : %d\n", int(psg->score_peel_sequence()));    
}

Peeler::Peeler(const Peeler& rhs) :
    ped(rhs.ped),
    map(rhs.map),
    rfunctions(rhs.rfunctions),
    trait_prob(rhs.trait_prob),
    lod_score(rhs.lod_score),
    initialised(rhs.initialised),
    locus(rhs.locus),
    count(rhs.count) {}

Peeler& Peeler::operator=(const Peeler& rhs) {
    
    if(&rhs != this) {
        ped = rhs.ped;
        map = rhs.map;
        rfunctions = rhs.rfunctions;
        trait_prob = rhs.trait_prob;
        lod_score = rhs.lod_score;
        initialised = rhs.initialised;
        locus = rhs.locus;
        count = rhs.count;
    }
    
    return *this;
}

Peeler::~Peeler() {}

double Peeler::get_trait_prob() {
    return trait_prob;
}

double Peeler::calc_trait_prob() {
    for(unsigned i = 0; i < rfunctions.size(); ++i) {
        rfunctions[i].evaluate(NULL, 0.5);
    }
    
    TraitRfunction& rf = rfunctions.back();
    
    return log(rf.get_result());
}

void Peeler::process(DescentGraph* dg) {

    for(unsigned i = 0; i < rfunctions.size(); ++i) {
        rfunctions[i].evaluate(dg, 0.5);
    }
    
    TraitRfunction& rf = rfunctions.back();
    double prob = log(rf.get_result()) - dg->get_recombination_prob(locus, false) - dg->get_marker_transmission();
    
    lod_score = initialised ? log_sum(prob, lod_score) : prob, initialised = true;
    ++count;
    
    if(locus == 0) {
    //printf("PROB %f\n", prob);
        fprintf(stderr, "%f %f\n", (lod_score - log(count) - trait_prob) / log(10.0), (prob - trait_prob) / log(10.0));
    }
}

double Peeler::get() {
    return (lod_score - log(count) - trait_prob) / log(10.0);
}

