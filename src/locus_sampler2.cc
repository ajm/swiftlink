using namespace std;

#include <vector>
#include <cmath>

#include "peel_sequence_generator.h"
#include "descent_graph_types.h"
#include "descent_graph.h"
#include "pedigree.h"
#include "genetic_map.h"
#include "rfunction_builder.h"
#include "sampler_rfunction.h"
#include "locus_sampler2.h"


void LocusSampler::init_rfunctions(PeelSequenceGenerator& psg) {
    //PeelSequenceGenerator psg(ped);
    //psg.build_peel_order();
    
    vector<PeelOperation>& ops = psg.get_peel_order();
    
    RfunctionBuilder<SamplerRfunction> build(ped, map, rfunctions);
    
    for(unsigned i = 0; i < ops.size(); ++i) {
        rfunctions.push_back(build.createRfunction(ops[i]));
    }
}

void LocusSampler::copy_rfunctions(const LocusSampler& rhs) {
    rfunctions.clear();
    
    for(unsigned i = 0; i < rhs.rfunctions.size(); ++i) {
        SamplerRfunction* rf = new SamplerRfunction(*(rhs.rfunctions[i]));
        rfunctions.push_back(rf);
    }
}

void LocusSampler::kill_rfunctions() {
    for(unsigned i = 0; i < rfunctions.size(); ++i) {
        delete rfunctions[i];
    }
}

unsigned LocusSampler::sample_hetero_mi(enum trait allele, enum phased_trait trait) {
    if(allele == TRAIT_U) {
        return (trait == TRAIT_UA) ? 0 : 1;
    }
    else {
        return (trait == TRAIT_UA) ? 1 : 0;
    }
}

// find prob of setting mi to 0
// find prob of setting mi to 1
// normalise + sample
unsigned LocusSampler::sample_homo_mi(DescentGraph& dg, unsigned personid, unsigned locus, enum parentage parent) {
    double prob_dist[2];
    double total;
    
    prob_dist[0] = 1.0;
    prob_dist[1] = 1.0;
    
    if(locus != 0) {
        prob_dist[0] *= ((dg.get(personid, locus - 1, parent) == 0) ? map->get_inversetheta(locus - 1) : map->get_theta(locus - 1));
        prob_dist[1] *= ((dg.get(personid, locus - 1, parent) == 1) ? map->get_inversetheta(locus - 1) : map->get_theta(locus - 1));
    }
    
    if(locus != (map->num_markers() - 1)) {
        prob_dist[0] *= ((dg.get(personid, locus + 1, parent) == 0) ? map->get_inversetheta(locus) : map->get_theta(locus));
        prob_dist[1] *= ((dg.get(personid, locus + 1, parent) == 1) ? map->get_inversetheta(locus) : map->get_theta(locus));
    }
    
    total = prob_dist[0] + prob_dist[1];
    
    prob_dist[0] /= total;
    prob_dist[1] /= total;
    
    
    return (get_random() < prob_dist[0]) ? 0 : 1;
}

// if a parent is heterozygous, then there is one choice of meiosis indicator
// if a parent is homozygous, then sample based on meiosis indicators to immediate left and right    
unsigned LocusSampler::sample_mi(DescentGraph& dg, \
                                 enum trait allele, enum phased_trait trait, \
                                 unsigned personid, unsigned locus, \
                                 enum parentage parent) {
    switch(trait) {
        case TRAIT_UA:
        case TRAIT_AU:
            return sample_hetero_mi(allele, trait);
            
        case TRAIT_UU:
        case TRAIT_AA:
            return sample_homo_mi(dg, personid, locus, parent);
            
        default:
            abort();
    }
}

// sample meiosis indicators
// if a parent is heterozygous, then there is one choice of meiosis indicator
// if a parent is homozygous, then sample based on meiosis indicators to immediate left and right
void LocusSampler::sample_meiosis_indicators(PeelMatrixKey& pmk, DescentGraph& dg, unsigned locus) {
    for(unsigned i = 0; i < ped->num_members(); ++i) {
        Person* p = ped->get_by_index(i);
        
        if(p->isfounder()) {
            continue;
        }
        
        enum phased_trait mat_trait = pmk.get(p->get_maternalid());
        enum phased_trait pat_trait = pmk.get(p->get_paternalid());
        
        enum phased_trait trait = pmk.get(i);
        enum trait mat_allele = ((trait == TRAIT_UU) or (trait == TRAIT_UA)) ? TRAIT_U : TRAIT_A;
        enum trait pat_allele = ((trait == TRAIT_UU) or (trait == TRAIT_AU)) ? TRAIT_U : TRAIT_A;
        
        unsigned mat_mi = sample_mi(dg, mat_allele, mat_trait, i, locus, MATERNAL);
        unsigned pat_mi = sample_mi(dg, pat_allele, pat_trait, i, locus, PATERNAL);
        
        dg.set(i, locus, MATERNAL, mat_mi);
        dg.set(i, locus, PATERNAL, pat_mi);        
    }
}

void LocusSampler::step(DescentGraph& dg, unsigned parameter) {
    unsigned locus = parameter;
    
    // forward peel
    for(unsigned i = 0; i < rfunctions.size(); ++i) {
        SamplerRfunction* rf = rfunctions[i];
        rf->evaluate(&dg, locus, 0.0);
    }
    
    PeelMatrixKey pmk;
    
    // reverse peel, sampling ordered genotypes
    for(int i = static_cast<int>(rfunctions.size()) - 1; i >= 0; --i) {
        SamplerRfunction* rf = rfunctions[i];
        rf->sample(pmk);
    }
    
    sample_meiosis_indicators(pmk, dg, locus);
    
    /*
    // XXX comment out when I know everything is cool
    if(not dg.likelihood()) {
        fprintf(stderr, "Error: descent graph produced by L-sampler is illegal!\n");
        abort();
    }
    */
}

void LocusSampler::reset() {
    set_all(false, false);
}

void LocusSampler::set_all(bool left, bool right) {
    for(unsigned i = 0; i < rfunctions.size(); ++i)
        rfunctions[i]->set_ignore(left, right);
}

void LocusSampler::sequential_imputation(DescentGraph& dg) {
    unsigned starting_locus = get_random_locus();
    
    set_all(true, true);
    step(dg, starting_locus);
    
    // iterate left through the markers
    set_all(true, false);
    for(int i = (starting_locus - 1); i >= 0; --i) {
        step(dg, i);
    }
    
    // iterate right through the markers
    set_all(false, true);
    for(int i = (starting_locus + 1); i < int(map->num_markers()); ++i) {
        step(dg, i);
    }
}
