using namespace std;

#include <cmath>
#include <vector>
#include <algorithm>

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
    //int lsampler_count = 0;
    //int msampler_count = 0;
    
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

    unsigned num_meioses = 2 * (ped->num_members() - ped->num_founders());
    
    vector<int> l_ordering;
    vector<int> m_ordering;
    
    for(unsigned int i = 0; i < map->num_markers(); ++i)
        l_ordering.push_back(i);
        
    for(unsigned int i = 0; i < num_meioses; ++i)
        m_ordering.push_back(i);
    
    random_shuffle(l_ordering.begin(), l_ordering.end());
    random_shuffle(m_ordering.begin(), m_ordering.end());
    /*
    for(unsigned int i = 0; i < num_meioses; ++i) {
        fprintf(stderr, "%d ", m_ordering[i]);
    }
    fprintf(stderr, "\n");
    */
    
    
    Progress p("MCMC: ", options.iterations + options.burnin);
    
    if(dg.get_likelihood() == LOG_ILLEGAL) {
        fprintf(stderr, "error: descent graph illegal pre-markov chain...\n");
        abort();
    }

        
    for(int i = 0; i < (options.iterations + options.burnin); ++i) {
        if(get_random() < options.lsampler_prob) {
            //lsampler_count ++;
        
            /*
            int j = get_random_int(map->num_markers());
            lsamplers[j]->step(dg, j);
            */
            
            
            for(unsigned int j = 0; j < map->num_markers(); ++j) {
                lsamplers[j]->step(dg, j);
            
                /*
                if(dg.get_likelihood() == LOG_ILLEGAL) {
                    fprintf(stderr, "error: descent graph illegal post-lsampler...\n");
                    abort();
                }
                */
            }
            
            
            /*
            int batches = 2;
            
            for(int j = 0; j < batches; ++j) {
                //#pragma omp parallel for
                for(int k = j; k < int(map->num_markers()); k += batches) {
                    lsamplers[k]->step(dg, k);
                }
            }
            */
            
            /*
            #pragma omp parallel for
            for(int j = 0; j < int(map->num_markers()); j += 2) {
                lsamplers[j]->step(dg, j);
            }
            
            #pragma omp parallel for
            for(int j = 1; j < int(map->num_markers()); j += 2) {
                lsamplers[j]->step(dg, j);
            }
            */
            
            
            //random_shuffle(l_ordering.begin(), l_ordering.end());
            /*
            for(unsigned int j = 0; j < l_ordering.size(); ++j) {
                lsamplers[l_ordering[j]]->step(dg, l_ordering[j]);
            }
            */
        }
        else {
            //msampler_count ++;
            
            
            //msampler.reset(dg);
            
            
            for(unsigned int j = 0; j < num_meioses; ++j) {
                msampler.step(dg, j);
                /*
                if(dg.get_likelihood() == LOG_ILLEGAL) {
                    fprintf(stderr, "error: descent graph illegal post-msampler...\n");
                    abort();
                }
                */
            }
            
            
            /*
            int j = get_random_int(num_meioses);
            //msampler.reset(dg);
            msampler.step(dg, j);
            */
            
            //random_shuffle(m_ordering.begin(), m_ordering.end());
            /*
            for(unsigned int j = 0; j < m_ordering.size(); ++j) {
                msampler.reset(dg);
                msampler.step(dg, m_ordering[j]);
            }
            */
        }
        
        p.increment();
        
        //if((i % options.scoring_period) == 0)
        //    fprintf(stderr, "%d %e\n", i, dg.get_likelihood() / log(10.0));
        
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
    
    /*
    printf("L-sampler : %.3f\nM-sampler : %.3f\n", \
        lsampler_count / double(lsampler_count + msampler_count), \
        msampler_count / double(lsampler_count + msampler_count));
    */
    
    return lod_scores;
}

