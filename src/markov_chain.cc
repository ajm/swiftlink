using namespace std;

#include <cmath>
#include <vector>
#include <algorithm>

#include <time.h>
#include <omp.h>

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
#include "random.h"


//#define MICROBENCHMARK_TIMING 1
//#define CODA_OUTPUT 1


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

    unsigned num_meioses = 2 * (ped->num_members() - ped->num_founders());
    
    vector<int> l_ordering;
    vector<int> m_ordering;
    
    /*
    for(unsigned int i = 0; i < map->num_markers(); ++i)
        l_ordering.push_back(i);
    */
    
    int markers_per_window = (map->num_markers() / omp_get_max_threads()) + \
                            ((map->num_markers() % omp_get_max_threads()) == 0 ? 0 : 1);
    fprintf(stderr, "setting %d markers per window\n", markers_per_window);
    //int markers_per_window = 98; // this needs to be set dynamically
    //int markers_per_window = 49;
    //int markers_per_window = 13;
    //int markers_per_window = 2;
    for(int i = 0; i < markers_per_window; ++i)
        l_ordering.push_back(i);
    
    
    for(unsigned int i = 0; i < num_meioses; ++i) {
        unsigned person_id = ped->num_founders() + (i / 2);
        enum parentage p = static_cast<enum parentage>(i % 2);
        
        Person* tmp = ped->get_by_index(person_id);
        
        if(not tmp->safe_to_ignore_meiosis(p)) {
            m_ordering.push_back(i);
        }
    }
    
    
    printf("P(l-sampler) = %f\n", options.lsampler_prob);
    
    Progress p("MCMC: ", options.iterations + options.burnin);
    
    if(dg.get_likelihood() == LOG_ILLEGAL) {
        fprintf(stderr, "error: descent graph illegal pre-markov chain...\n");
        abort();
    }

        
    for(int i = 0; i < (options.iterations + options.burnin); ++i) {
        if(get_random() < options.lsampler_prob) {
            /*
            random_shuffle(l_ordering.begin(), l_ordering.end());
            
            for(unsigned int j = 0; j < l_ordering.size(); ++j) {
                lsamplers[l_ordering[j]]->step(dg, l_ordering[j]);
            }
            */
            
            random_shuffle(l_ordering.begin(), l_ordering.end());
            
            // run every nth locus indepentently, shuffling the order of n
            for(unsigned int j = 0; j < l_ordering.size(); ++j) {
                #pragma omp parallel for
                for(unsigned int k = l_ordering[j]; k < map->num_markers(); k += markers_per_window) {
                    lsamplers[k]->step(dg, k);
                }
            }
        }
        else {
            random_shuffle(m_ordering.begin(), m_ordering.end());
            
            msampler.reset(dg, m_ordering[0]);
            for(unsigned int j = 0; j < m_ordering.size(); ++j) {
                msampler.step(dg, m_ordering[j]);
            }
        }
        
        
        p.increment();
        
        /*
        double current_likelihood = dg.get_likelihood();
        if(current_likelihood == LOG_ILLEGAL) {
            fprintf(stderr, "error: descent graph illegal...\n");
            abort();
        }
        */
        
        #ifdef CODA_OUTPUT
        if(i < options.burnin) {
            double current_likelihood = dg.get_likelihood();
            if(current_likelihood == LOG_ILLEGAL) {
                fprintf(stderr, "error: descent graph illegal...\n");
                abort();
            }
            
            fprintf(stderr, "%d\t%f\n", i+1, current_likelihood);
        }
        #endif
        
        if(i < options.burnin) {
            continue;
        }
                
        if((i % options.scoring_period) == 0) {
        
            #ifdef MICROBENCHMARK_TIMING
            struct timespec start_time;
            struct timespec end_time;
            
            clock_gettime(CLOCK_MONOTONIC, &start_time);
            
            int repeats = 100;
            
            for(int x = 0; x < repeats; ++x) {
            #endif
            
            #pragma omp parallel for
            for(int j = 0; j < int(map->num_markers() - 1); ++j) {
                peelers[j]->process(&dg);
            }
            
            #ifdef MICROBENCHMARK_TIMING
            }
            
            clock_gettime(CLOCK_MONOTONIC, &end_time);
            
            double milliseconds = \
                (((end_time.tv_sec   * 1000.0) + (end_time.tv_nsec   / 1000000.0)) - \
                 ((start_time.tv_sec * 1000.0) + (start_time.tv_nsec / 1000000.0))) / double(repeats);
            
            fprintf(stderr, "LODSCORE %.3f\n", milliseconds);
            #endif
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

