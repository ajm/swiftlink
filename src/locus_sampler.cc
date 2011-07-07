using namespace std;

#include <vector>
#include <cmath>

#include "descent_graph_types.h"
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
#include "progress.h"
#include "peeler.h"

#include "linkage_writer.h"

/*
    this is pretty confusing, perhaps in need of a refactor(?)
    
    TraitRfunction deals with disease traits
    SamplerRfunction deals with phased genotypes
        both disease traits and phased genotypes are represented by "enum phased_trait"
    DescentGraph deals with meiosis indicators
        these are just zeros and ones
        
 */

LocusSampler::LocusSampler(Pedigree* ped, GeneticMap* map) :
    ped(ped), 
    map(map), 
    dg(ped, map),
    rfunctions(),
    peel(ped, map),
    burnin_steps(0) {
    
    dg.random_descentgraph();
    
    init_rfunctions();
}

LocusSampler::LocusSampler(const LocusSampler& rhs) :
    ped(ped),
    map(map),
    dg(dg),
    rfunctions(),
    peel(rhs.peel),
    burnin_steps(rhs.burnin_steps) {
    
    copy_rfunctions(rhs);
}

LocusSampler& LocusSampler::operator=(const LocusSampler& rhs) {
    if(this != &rhs) {
		ped = rhs.ped;
		map = rhs.map;
		dg = rhs.dg;
		peel = rhs.peel;
		burnin_steps = rhs.burnin_steps;
		
		kill_rfunctions();
		copy_rfunctions(rhs);
	}
	        
	return *this;
}

LocusSampler::~LocusSampler() {
    kill_rfunctions();
}

void LocusSampler::init_rfunctions() {
    PeelSequenceGenerator psg(ped);
    psg.build_peel_order();
    
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

void LocusSampler::anneal(unsigned iterations) {
    unsigned steps_per_temp = iterations / 1000;
    double annealing_temp = 1.0;
    
    Progress q("Simulated Annealing:", iterations);
    q.start();
    
    for(unsigned i = 0; i < iterations; ++i) {
        step(annealing_temp);
        
        q.increment();
        
        if((i % steps_per_temp) == 0) {
            annealing_temp *= 0.99;
        }
    }
    
    q.finish();
}

double LocusSampler::get_random() {
    return random() / static_cast<double>(RAND_MAX);
}

unsigned LocusSampler::get_random(unsigned i) {
    return random() % i;
}

unsigned LocusSampler::get_random_locus() {
    return get_random(ped->num_markers());
}

unsigned LocusSampler::sample_hetero_mi(unsigned allele, enum phased_trait trait) {
    if(allele == 0) {
        return (trait == TRAIT_UA) ? 0 : 1;
    }
    else {
        return (trait == TRAIT_UA) ? 1 : 0;
    }
}

unsigned LocusSampler::sample_homo_mi(unsigned personid, unsigned locus, enum parentage parent, double temperature) {
    // find prob of setting mi to 0
    // find prob of setting mi to 1
    // normalise + sample
    double prob_dist[2];
    double theta;
    double total;
    
    prob_dist[0] = 1.0;
    prob_dist[1] = 1.0;
    
    if(locus != 0) {
        theta = exp(map->get_theta(locus - 1, temperature));
        prob_dist[0] *= ((dg.get(personid, locus - 1, parent) == 0) ? 1.0 - theta : theta);
        prob_dist[1] *= ((dg.get(personid, locus - 1, parent) == 1) ? 1.0 - theta : theta);
    }
    
    if(locus != (map->num_markers() - 1)) {
        theta = exp(map->get_theta(locus, temperature));
        prob_dist[0] *= ((dg.get(personid, locus + 1, parent) == 0) ? 1.0 - theta : theta);
        prob_dist[1] *= ((dg.get(personid, locus + 1, parent) == 1) ? 1.0 - theta : theta);
    }
    
    total = prob_dist[0] + prob_dist[1];
    
    prob_dist[0] /= total;
    prob_dist[1] /= total;
    
    
    return (get_random() < prob_dist[0]) ? 0 : 1;
}

// if a parent is heterozygous, then there is one choice of meiosis indicator
// if a parent is homozygous, then sample based on meiosis indicators to immediate left and right    
unsigned LocusSampler::sample_mi(unsigned allele, enum phased_trait trait, \
                                 unsigned personid, unsigned locus, enum parentage parent, double temperature) {
    switch(trait) {
        case TRAIT_UA:
        case TRAIT_AU:
            return sample_hetero_mi(allele, trait);
        
        case TRAIT_UU:
        case TRAIT_AA:
            return sample_homo_mi(personid, locus, parent, temperature);
        
        default:
            abort();
    }
}

void LocusSampler::step(double temperature) {
    unsigned locus = get_random_locus();
    
    // forward peel
    for(unsigned i = 0; i < rfunctions.size(); ++i) {
        SamplerRfunction* rf = rfunctions[i];
        rf->evaluate(&dg, locus, 0.0, temperature);
    }
    
    PeelMatrixKey pmk;
    
    // reverse peel
    for(int i = static_cast<int>(rfunctions.size()) - 1; i > -1; --i) {
        SamplerRfunction* rf = rfunctions[i];
        rf->sample(pmk);
    }
    
    // sample meiosis indicators
    // if a parent is heterozygous, then there is one choice of meiosis indicator
    // if a parent is homozygous, then sample based on meiosis indicators to immediate left and right
    for(unsigned i = 0; i < ped->num_members(); ++i) {
        Person* p = ped->get_by_index(i);
        
        if(p->isfounder()) {
            continue;
        }
        
        enum phased_trait mat_trait = pmk.get(p->get_maternalid());
        enum phased_trait pat_trait = pmk.get(p->get_paternalid());
        
        enum phased_trait trait = pmk.get(i);
        unsigned mat_allele = ((trait == TRAIT_UU) or (trait == TRAIT_UA)) ? 0 : 1;
        unsigned pat_allele = ((trait == TRAIT_UU) or (trait == TRAIT_AU)) ? 0 : 1;
        
        unsigned mat_mi = sample_mi(mat_allele, mat_trait, i, locus, MATERNAL, temperature);
        unsigned pat_mi = sample_mi(pat_allele, pat_trait, i, locus, PATERNAL, temperature);
        
        // just do a set on the DescentGraph?
        dg.set(i, locus, MATERNAL, mat_mi);
        dg.set(i, locus, PATERNAL, pat_mi);
    }
}

void LocusSampler::run(unsigned start_step, unsigned iterations, double temperature, Peeler& p) {
    unsigned end_step = start_step + iterations;
    
    for(unsigned i = start_step; i < end_step; ++i) {
        step(temperature);
        
        if((temperature != 0.0) or (i < burnin_steps)) {
            continue;
        }
        
        if((i % SAMPLING_PERIOD) == 0) {
            p.process(dg);
        }
    }
}

unsigned LocusSampler::update_temperature_hastings(unsigned temps, unsigned current_temp) {
    unsigned new_temp;
    double old_prob;
    double new_prob;
    double prob = 1.0;
    
    if(current_temp == 0) {
        new_temp = 1;
    }
    else if(current_temp == temps) {
        new_temp = temps - 1;
    }
    else {
        new_temp = current_temp + ((get_random() < 0.5) ? 1 : -1); 
        prob = 0.5;
    }
    
    dg.likelihood(&old_prob, current_temp / static_cast<double>(temps+1));
    dg.likelihood(&new_prob, new_temp     / static_cast<double>(temps+1));
    
    double hastings_ratio = log(0.5 / prob);
    
    return (log(get_random()) < ((new_prob - old_prob) + hastings_ratio)) ? \
        new_temp : current_temp;
}

unsigned LocusSampler::update_temperature(unsigned temps, unsigned current_temp) {
    
    int diff = get_random() < 0.5 ? 1 : -1 ;
    int new_temp = current_temp + diff;
    
    //fprintf(stderr, "Z %d\n", diff);
    
    if(new_temp < 0) {
        new_temp = 0;
    }
    
    if(new_temp > static_cast<int>(temps)) {
        new_temp = temps;
    }
    
    double old_prob;
    double new_prob;
    
    dg.likelihood(&old_prob, current_temp / static_cast<double>(temps+1));
    dg.likelihood(&new_prob, new_temp     / static_cast<double>(temps+1));
    
    return (log(get_random()) < (new_prob - old_prob)) ? static_cast<unsigned>(new_temp) : current_temp;
}

Peeler* LocusSampler::temper(unsigned iterations, unsigned temperatures) {
    unsigned burnin_steps = iterations * 0.1;
    unsigned temperature_max = temperatures - 1;
    unsigned temperature_level = temperature_max;
    double theta = 0.0;
    unsigned num_samples = 0;
    
    unsigned* temp_counts = new unsigned[temperatures];
    for(unsigned i = 0; i < temperatures; ++i)
        temp_counts[i] = 0;
    
    
    unsigned steps_per_temp = iterations / 800;
    Progress q("Simulated Annealing:", iterations);
    q.start();
    theta = 1.0;
    temperature_level = temperature_max;
    for(unsigned i = 0; i < iterations; ++i) {
        step(theta);
        
        q.increment();
        
        if((i % steps_per_temp) == 0) {
            theta *= 0.99;
        }
    }
    
    q.finish();
    
    
    theta = 0.0;
    temperature_level = 0;
    //theta = 1.0;
    
    //iterations *= 10;
    //burnin_steps = iterations * 0.1;
    
    Progress p("Simulated Tempering:", iterations);
    p.start();
    
    
    for(unsigned i = 0; i < iterations; ++i) {
        step(theta);
        p.increment();
        
        if((i % 5) == 0) {
            temperature_level = update_temperature_hastings(temperature_max, temperature_level);
            theta = temperature_level / static_cast<double>(temperatures);
            
            temp_counts[temperature_level] += 1;
        }
        
        if(i < burnin_steps) {    
            continue;
        }
        
        if((i % 5) == 0) {
            if(temperature_level == 0) {
                peel.process(dg);
                num_samples++;
                
                printf("\n\n\n%d samples\n", num_samples);
                LinkageWriter lw(map, &peel, "x", true);
                lw.write();
            }
        }
    }
    
    p.finish();
    
    //fprintf(stdout, "total samples = %d, count at temperature 0 = %d\n", num_samples, num_temp0);
    fprintf(stderr, "#samples = %d\n", num_samples);
    
    
    for(unsigned i = 0; i < temperatures; ++i)
        fprintf(stderr, "%d %d\n", i, temp_counts[i]);
    
    delete[] temp_counts;
    
    return &peel;
}

double LocusSampler::likelihood(double temperature) {
    double tmp;
    
    dg.likelihood(&tmp, temperature);
    
    return tmp;
}

