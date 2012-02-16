using namespace std;

#include <cstdio>

#include "sequential_imputation.h"
#include "logarithms.h"
#include "progress.h"
#include "descent_graph.h"
#include "locus_sampler2.h"
#include "peel_sequence_generator.h"

/*
void SequentialImputation::run(DescentGraph& dg, int iterations) {
    DescentGraph tmp(ped, map);
    double tmp_prob, best_prob = LOG_ILLEGAL;
    
    LocusSampler lsampler(ped, map, psg, 0);
    
    tmp.random_descentgraph();
    
    if((tmp_prob = tmp.get_likelihood()) == LOG_ZERO) {
        fprintf(stderr, "error: failed to produce a valid random descent graph\n");
        abort();
    }
    
    //fprintf(stderr, "initial random likelihood = %e\n", tmp_prob);
    
    Progress p("Sequential Imputation: ", iterations);
    
    do {
        lsampler.sequential_imputation(tmp);
    
        if((tmp_prob = tmp.get_likelihood()) == LOG_ZERO) {
            p.finish();
            fprintf(stderr, "error: sequential imputation produced an invalid descent graph\n");
            fprintf(stderr, "%s\n", tmp.debug_string().c_str());
            abort();
        }
        
        p.increment();
        
        if(tmp_prob > best_prob) {
            dg = tmp;
            best_prob = tmp_prob;
        }
        
    } while(--iterations > 0);
    
    p.finish();
    
    //fprintf(stderr, "starting point likelihood = %e\n", best_prob);
}
*/

void SequentialImputation::run(DescentGraph& dg, int iterations) {
    double best_prob = LOG_ZERO;
    
    if(iterations == 0) {
        return;
    }
    
    vector<LocusSampler*> lsamplers;
    vector<DescentGraph*> graphs;
    vector<double> likelihoods(iterations, 0.0);
    for(int i = 0; i < iterations; ++i) {
        LocusSampler* tmp = new LocusSampler(ped, map, psg, 0);
        lsamplers.push_back(tmp);
        
        DescentGraph* tmp_dg = new DescentGraph(ped, map);
        graphs.push_back(tmp_dg);
    }
    
    Progress p("Sequential Imputation: ", iterations);
    
    #pragma omp parallel for
    for(int i = 0; i < iterations; ++i) {
        likelihoods[i] = lsamplers[i]->sequential_imputation(*(graphs[i]));
        /*
        if(graphs[i]->get_likelihood() == LOG_ZERO) {
            fprintf(stderr, "error: sequential imputation produced an invalid descent graph\n");
            fprintf(stderr, "%s\n", graphs[i]->debug_string().c_str());
            abort();
        }
        */
        p.increment();
    }
    
    p.finish();
    
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
        
        //printf("likelihood = %e (%e)\n", likelihoods[i], likelihoods[i] / log(10));
    }
    
    dg = *graphs[index];
    
    printf("starting point likelihood = %e (%e)\n", best_prob, best_prob / log(10));
    
    for(int i = 0; i < iterations; ++i) {
        delete lsamplers[i];
        delete graphs[i];
    }
}

