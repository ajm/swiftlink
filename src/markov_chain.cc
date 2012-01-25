using namespace std;

#include <cmath>

#include "markov_chain.h"
#include "peel_sequence_generator.h"
#include "peeler.h"
#include "locus_sampler2.h"
#include "meiosis_sampler.h"
#include "pedigree.h"
#include "genetic_map.h"
#include "descent_graph.h"
#include "linkage_writer.h"
#include "progress.h"
#include "types.h"


double* MarkovChain::run(DescentGraph& dg) {
    // lod scorers
    vector<Peeler*> peelers;
    for(unsigned int i = 0; i < (map->num_markers() - 1); ++i) {
        Peeler* tmp = new Peeler(ped, map, psg, i);
        peelers.push_back(tmp);
    }
    
    // create samplers
    vector<LocusSampler*> lsamplers;
    for(unsigned int i = 0; i < map->num_markers(); ++i) {
        LocusSampler* tmp = new LocusSampler(ped, map, psg, i);
        lsamplers.push_back(tmp);
    }
    
    MeiosisSampler msampler(ped, map);
    
    
    Progress p("MCMC: ", options.iterations + options.burnin);
    
    
    unsigned num_meioses = 2 * (ped->num_members() - ped->num_founders());
        
    for(int i = 0; i < (options.iterations + options.burnin); ++i) {
        if(get_random() < options.lsampler_prob) {
            #pragma omp parallel for
            for(int j = 0; j < int(map->num_markers()); j += 2) {
                lsamplers[j]->step(dg, j);
            }
            
            #pragma omp parallel for
            for(int j = 1; j < int(map->num_markers()); j += 2) {
                lsamplers[j]->step(dg, j);
            }
            
            /*
            for(unsigned int j = 0; j < map->num_markers(); ++j) {
                lsamplers[j]->step(dg, j);
            }
            */
        }
        else {
            msampler.reset(dg);
            for(unsigned int j = 0; j < num_meioses; ++j) {
                msampler.step(dg, j);
            }
        }
        
        p.increment();
        
        if(i < options.burnin) {
            continue;
        }
        
        if((i % options.scoring_period) == 0) {
            #pragma omp parallel for
            for(int j = 0; j < int(map->num_markers() - 1); ++j) {
                peelers[j]->process(&dg);
            }
        }
    }
    
    p.finish();
    
    
    // dump lsamplers
    for(unsigned i = 0; i < lsamplers.size(); ++i) {
        delete lsamplers[i];
    }
    
    double* lod_scores = new double[map->num_markers() - 1];
    for(unsigned int i = 0; i < (map->num_markers() - 1); ++i) {
        lod_scores[i] = peelers[i]->get();
        delete peelers[i];
    }
    
    return lod_scores;
}

