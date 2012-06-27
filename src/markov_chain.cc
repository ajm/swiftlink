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
#include "lod_score.h"
#include "gpu_lodscores.h"


//#define MICROBENCHMARK_TIMING 1
//#define CODA_OUTPUT 1


LODscores* MarkovChain::run(DescentGraph& dg) {

    LODscores* lod = new LODscores(map);
    GPULodscores* gpulod = 0;
    
    // lod scorers
    vector<Peeler*> peelers;
    for(int i = 0; i < omp_get_max_threads(); ++i) {
        Peeler* tmp = new Peeler(ped, map, psg, lod);
        //tmp->set_locus(i);
        peelers.push_back(tmp);
    }
    
    double trait_prob = peelers[0]->calc_trait_prob();
    
    lod->set_trait_prob(trait_prob);
    
    printf("P(T) = %e\n", trait_prob / log(10));
    
    exit(-1);
    
    if(options.use_gpu) {
        gpulod = new GPULodscores(ped, map, psg, options, trait_prob);
    }
    
    
    // create samplers
    vector<LocusSampler*> lsamplers;
    for(int i = 0; i < omp_get_max_threads(); ++i) {
        LocusSampler* tmp = new LocusSampler(ped, map, psg, i);
        lsamplers.push_back(tmp);
    }
    
    
    
    MeiosisSampler msampler(ped, map);

    unsigned num_meioses = 2 * (ped->num_members() - ped->num_founders());
    
    vector<int> l_ordering;
    vector<int> m_ordering;
    
    
    int markers_per_window = (map->num_markers() / omp_get_max_threads()) + \
                            ((map->num_markers() % omp_get_max_threads()) == 0 ? 0 : 1);
    fprintf(stderr, "setting %d markers per window\n", markers_per_window);
    
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
            
            random_shuffle(l_ordering.begin(), l_ordering.end());
            
            // run every nth locus indepentently, shuffling the order of n
            for(unsigned int j = 0; j < l_ordering.size(); ++j) {
                #pragma omp parallel for
                for(unsigned int k = l_ordering[j]; k < map->num_markers(); k += markers_per_window) {
                    lsamplers[omp_get_thread_num()]->set_locus_minimal(k);
                    lsamplers[omp_get_thread_num()]->step(dg, k);
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
         
            if(not options.use_gpu) {
                #pragma omp parallel for
                for(int j = 0; j < int(map->num_markers() - 1); ++j) {
                    peelers[omp_get_thread_num()]->set_locus(j);
                    peelers[omp_get_thread_num()]->process(&dg);
                }
            }
            else {
                gpulod->calculate(dg);
            }
        }
        
    }
    
    
    p.finish();
    
    if(options.use_gpu) {
        gpulod->get_results(lod);
    }
    
    // dump lsamplers
    for(unsigned int i = 0; i < lsamplers.size(); ++i) {
        delete lsamplers[i];
    }
    
    for(unsigned int i = 0; i < peelers.size(); ++i) {
        delete peelers[i];
    }
    
    return lod;
}

