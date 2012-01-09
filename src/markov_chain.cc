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
//#include "gpu_wrapper.h"


void MarkovChain::initialise(DescentGraph& dg, PeelSequenceGenerator& psg) {
    DescentGraph tmp(ped, map);
    double tmp_prob, best_prob = LOG_ILLEGAL;
    int iterations = 100;
    //int iterations = 1;
    
    LocusSampler lsampler(ped, map, psg, 0);
    
    tmp.random_descentgraph();
    if((tmp_prob = tmp.get_likelihood()) == LOG_ZERO) {
        fprintf(stderr, "error: failed to produce a valid random descent graph\n");
        abort();
    }
    
    fprintf(stderr, "initial random likelihood = %e\n", tmp_prob);
    
    do {
        lsampler.sequential_imputation(tmp);
    
        if((tmp_prob = tmp.get_likelihood()) == LOG_ZERO) {
            fprintf(stderr, "error: sequential imputation produced an invalid descent graph\n");
            fprintf(stderr, "%s\n", tmp.debug_string().c_str());
            abort();
        }
        
        //fprintf(stderr, "sequential imputation = %e\n", tmp_prob);
        
        if(tmp_prob > best_prob) {
            dg = tmp;
            best_prob = tmp_prob;
        }
        
    } while(--iterations > 0);
    
    fprintf(stderr, "starting point likelihood = %e\n", best_prob);
}

Peeler* MarkovChain::run(unsigned iterations, double temperature) {
    unsigned burnin = iterations * 0.1;
    
    //iterations = 1020;
    //burnin = 1000;
    
    map->set_temperature(temperature);

    // create a descent graph
    DescentGraph dg(ped, map);
    
    // build peeling sequence for L-sampler and Peeler
    PeelSequenceGenerator psg(ped);
    psg.build_peel_order();
    
    initialise(dg, psg);
    
    // create an object to perform LOD scoring
    // allocate on the heap so we can return it
    Peeler* peel = new Peeler(ped, map, psg);
    
    // create samplers
    vector<LocusSampler*> lsamplers;
    for(unsigned int i = 0; i < map->num_markers(); ++i) {
        LocusSampler* tmp = new LocusSampler(ped, map, psg, i);
        lsamplers.push_back(tmp);
    }
    
    MeiosisSampler msampler(ped, map);
    
    /*
    GPUWrapper gpu(ped, map, psg);
    
    gpu.run(dg, iterations, burnin, 10, peel->get_trait_prob());
    
    exit(0);
    */
    //exit(0);
    
    Progress p("MCMC: ", iterations);
    
    
    unsigned num_meioses = 2 * (ped->num_members() - ped->num_founders());
        
    for(unsigned int i = 0; i < iterations; ++i) {
        if((random() / DBL_RAND_MAX) < 0.5) {
            /*
            for(unsigned int j = 0; j < map->num_markers(); j += 2) {
                lsamplers[j]->step(dg, j);
            }
            for(unsigned int j = 1; j < map->num_markers(); j += 2) {
                lsamplers[j]->step(dg, j);
            }
            */
            for(unsigned int j = 0; j < map->num_markers(); ++j) {
                lsamplers[j]->step(dg, j);
            }
        }
        else {
            msampler.reset(dg);
            for(unsigned int j = 0; j < num_meioses; ++j) {
                msampler.step(dg, j);
            }
        }
        
        p.increment();
        
        if(i < burnin) {
            continue;
        }
        
        if((i % 10) == 0) {
            peel->process(dg);
        }
    }
    
    p.finish();
    
    
    // dump lsamplers
    for(unsigned i = 0; i < lsamplers.size(); ++i) {
        delete lsamplers[i];
    }
    
    return peel;
}
