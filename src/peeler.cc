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
#include "lod_score.h"


Peeler::Peeler(Pedigree* p, GeneticMap* g, PeelSequenceGenerator* psg, LODscores* lod) : 
    ped(p), 
    map(g),
    lod(lod),
    rfunctions(),
    locus(0) {
    
    vector<PeelOperation>& ops = psg->get_peel_order();
    
    rfunctions.reserve(ops.size()); // need to do this otherwise pointers may not work later...
    
    for(unsigned i = 0; i < ops.size(); ++i) {
        Rfunction* prev1 = ops[i].get_previous_op1() == -1 ? NULL : &rfunctions[ops[i].get_previous_op1()];
        Rfunction* prev2 = ops[i].get_previous_op2() == -1 ? NULL : &rfunctions[ops[i].get_previous_op2()];
        
        rfunctions.push_back(TraitRfunction(ped, map, locus, &(ops[i]), prev1, prev2));
    }   
}

Peeler::Peeler(const Peeler& rhs) :
    ped(rhs.ped),
    map(rhs.map),
    lod(rhs.lod),
    rfunctions(rhs.rfunctions),
    locus(rhs.locus) {}

Peeler& Peeler::operator=(const Peeler& rhs) {
    
    if(&rhs != this) {
        ped = rhs.ped;
        map = rhs.map;
        lod = rhs.lod;
        rfunctions = rhs.rfunctions;
        locus = rhs.locus;
    }
    
    return *this;
}

Peeler::~Peeler() {}

double Peeler::calc_trait_prob() {
    for(unsigned i = 0; i < rfunctions.size(); ++i) {
        rfunctions[i].evaluate(NULL, 0.5);
    }
    
    TraitRfunction& rf = rfunctions.back();
    
    return log(rf.get_result());
}

double Peeler::get_trait_prob() {
    return calc_trait_prob();
}

void Peeler::process(DescentGraph* dg) {

    for(unsigned i = 0; i < rfunctions.size(); ++i) {
        rfunctions[i].evaluate(dg, 0.5);
    }
    
    TraitRfunction& rf = rfunctions.back();
    
    double prob = log(rf.get_result()) - dg->get_recombination_prob(locus, false) - dg->get_marker_transmission();
     
    lod->add(locus, 0, prob);
}

