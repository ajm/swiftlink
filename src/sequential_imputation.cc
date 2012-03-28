using namespace std;

#include <cstdio>

#include "sequential_imputation.h"
#include "logarithms.h"
#include "progress.h"
#include "descent_graph.h"
#include "locus_sampler2.h"
#include "peel_sequence_generator.h"


void SequentialImputation::run(DescentGraph& dg, int iterations) {
    DescentGraph tmp(ped, map);
    double tmp_prob, best_prob = LOG_ZERO;
    
    LocusSampler lsampler(ped, map, psg, 0);
    
    Progress p("Sequential Imputation: ", iterations);
    
    do {
        tmp_prob = lsampler.sequential_imputation(tmp);
        
        if(tmp_prob == LOG_ZERO) {
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
    
    fprintf(stderr, "starting likelihood (log10) = %.3f\n", best_prob / log(10));
}


void SequentialImputation::parallel_run(DescentGraph& dg, int iterations) {
    double best_prob = LOG_ZERO;
    
    if(iterations == 0) {
        
        LocusSampler tmp(ped, map, psg, 0);
        double tmp_likelihood = tmp.locus_by_locus(dg);
        
        printf("starting likelihood (log10) = %.3f\n", tmp_likelihood / log(10));
        
        return;
    }
    
    vector<LocusSampler*> lsamplers;
    vector<DescentGraph*> graphs;
    vector<double> likelihoods(64, 0.0);
    for(int i = 0; i < 64; ++i) {
        LocusSampler* tmp = new LocusSampler(ped, map, psg, 0);
        lsamplers.push_back(tmp);
        
        DescentGraph* tmp_dg = new DescentGraph(ped, map);
        graphs.push_back(tmp_dg);
    }
    
    Progress p("Sequential Imputation: ", iterations);
    
    for(int i = 0; i < iterations; i += 64) {
        
        #pragma omp parallel for
        for(int j = 0; j < 64; ++j) {
            if((i + j) < iterations) {
                likelihoods[j] = lsamplers[j]->sequential_imputation(*(graphs[j]));
                p.increment();
            }
        }
        
        int index = -1;
        
        for(int j = 0; j < 64; ++j) {
            if((i + j) < iterations) {
                if(likelihoods[j] == LOG_ZERO) {
                    fprintf(stderr, "error: sequential imputation produced an invalid descent graph\n");
                    abort();
                }
                
                if(likelihoods[j] > best_prob) {
                    best_prob = likelihoods[j];
                    index = j;
                }
            }
        }
        
        if(index != -1) {
            dg = *graphs[index];
        }
    }
    
    p.finish();
    
    
    printf("starting likelihood (log10) = %.3f\n", best_prob / log(10));
    
    for(int i = 0; i < 64; ++i) {
        delete lsamplers[i];
        delete graphs[i];
    }
}

