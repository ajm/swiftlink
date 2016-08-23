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

using namespace std;


Peeler::Peeler(Pedigree* p, GeneticMap* g, PeelSequenceGenerator* psg, LODscores* lod, bool sex_linked) :
    ped(p), 
    map(g),
    lod(lod),
    rfunctions(),
    locus(0),
    sex_linked(sex_linked) {
    
    vector<PeelOperation>& ops = psg->get_peel_order();
    
    rfunctions.reserve(ops.size()); // need to do this otherwise pointers may not work later...
    
    for(unsigned int i = 0; i < ops.size(); ++i) {
        vector<unsigned int>& prev_indices = ops[i].get_prevfunctions();
        vector<Rfunction*> prev_pointers;
        
        for(unsigned int j = 0; j < prev_indices.size(); ++j) {
            prev_pointers.push_back(&(rfunctions[prev_indices[j]]));
        }
        
        rfunctions.push_back(TraitRfunction(ped, map, locus, &(ops[i]), prev_pointers, sex_linked));
    }
}

Peeler::Peeler(const Peeler& rhs) :
    ped(rhs.ped),
    map(rhs.map),
    lod(rhs.lod),
    rfunctions(rhs.rfunctions),
    locus(rhs.locus),
    sex_linked(rhs.sex_linked) {}

Peeler& Peeler::operator=(const Peeler& rhs) {
    
    if(&rhs != this) {
        ped = rhs.ped;
        map = rhs.map;
        lod = rhs.lod;
        rfunctions = rhs.rfunctions;
        locus = rhs.locus;
        sex_linked = rhs.sex_linked;
    }
    
    return *this;
}

Peeler::~Peeler() {}

double Peeler::calc_trait_prob() {
    for(unsigned i = 0; i < rfunctions.size(); ++i) {
        rfunctions[i].evaluate(NULL, 1);
    }
    
    TraitRfunction& rf = rfunctions.back();
    
    return log(rf.get_result());
}

double Peeler::get_trait_prob() {
    return calc_trait_prob();
}

void Peeler::process(DescentGraph* dg) {

    unsigned int num_lod_scores = lod->get_lodscores_per_marker();
    
    for(unsigned int i = 0; i < num_lod_scores; ++i) {
    
        for(unsigned j = 0; j < rfunctions.size(); ++j) {
            rfunctions[j].set_thetas(i + 1); // XXX really messy, change api somehow?
            rfunctions[j].evaluate(dg, i + 1); // geneticmap expects indexing from 1
        }
        
        TraitRfunction& rf = rfunctions.back();
        
        if(rf.get_result() <= 0.0) {
            fprintf(stderr, "\n\nerror: intermediate state had a likelihood of 0.0 or less (lod score %d, likelihood = %e)\n", i, rf.get_result());
            exit(1);
        }

        double prob = log(rf.get_result()) - \
                      dg->get_recombination_prob(locus, false) - \
                      dg->get_marker_transmission();
        
        lod->add(locus, i, prob);
    }
}

