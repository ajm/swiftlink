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
#include "gpu_wrapper.h"


void MarkovChain::initialise(DescentGraph& dg, PeelSequenceGenerator& psg) {
    DescentGraph tmp(ped, map);
    double tmp_prob, best_prob = LOG_ILLEGAL;
    //int iterations = 100;
    int iterations = 1;
    
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

void MarkovChain::parallel_initialise(DescentGraph& dg, PeelSequenceGenerator& psg) {
    double best_prob = LOG_ILLEGAL;
    int iterations = 100;
    
    vector<LocusSampler*> lsamplers;
    vector<DescentGraph*> graphs;
    vector<double> likelihoods(iterations, 0.0);
    for(int i = 0; i < iterations; ++i) {
        LocusSampler* tmp = new LocusSampler(ped, map, psg, 0);
        lsamplers.push_back(tmp);
        
        DescentGraph* tmp_dg = new DescentGraph(ped, map);
        graphs.push_back(tmp_dg);
    }
    
    #pragma omp parallel for
    for(int i = 0; i < iterations; ++i) {
        lsamplers[i]->sequential_imputation(*(graphs[i]));
        likelihoods[i] = graphs[i]->get_likelihood();
    }
    
    int index = 0;
    for(int i = 0; i < iterations; ++i) {
        if(likelihoods[i] == LOG_ZERO) {
            fprintf(stderr, "error: sequential imputation produced an invalid descent graph\n");
            fprintf(stderr, "%s\n", graphs[i]->debug_string().c_str());
            abort();
        }
        
        if(likelihoods[i] > best_prob) {
            best_prob = likelihoods[i];
            index = i;
        }
    }
    
    dg = *graphs[index];
    
    fprintf(stderr, "starting point likelihood = %e\n", best_prob);
    
    for(int i = 0; i < iterations; ++i) {
        delete lsamplers[i];
        delete graphs[i];
    }
}

double* MarkovChain::run(unsigned iterations, double temperature) {
    unsigned burnin = iterations * 0.1;
    
    //iterations = 1020;
    //burnin = 1000;
    
    map->set_temperature(temperature);

    // create a descent graph
    DescentGraph dg(ped, map);
    
    // build peeling sequence for L-sampler and Peeler
    PeelSequenceGenerator psg(ped);
    psg.build_peel_order();
    
    //initialise(dg, psg);
    parallel_initialise(dg, psg);
    
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
    
    
    GPUWrapper gpu(ped, map, psg);
    
    gpu.run(dg, iterations, burnin, 10, peelers[0]->get_trait_prob());
    
    exit(0);
    
    //exit(0);
    
    Progress p("MCMC: ", iterations);
    
    
    unsigned num_meioses = 2 * (ped->num_members() - ped->num_founders());
        
    for(unsigned int i = 0; i < iterations; ++i) {
        if((random() / DBL_RAND_MAX) < 0.5) {
            
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
        
        if(i < burnin) {
            continue;
        }
        
        if((i % 10) == 0) {
            #pragma omp parallel for
            for(int i = 0; i < int(map->num_markers() - 1); ++i) {
                peelers[i]->process(&dg);
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
