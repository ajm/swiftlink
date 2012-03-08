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
#include "random.h"


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
    for(unsigned int i = 0; i < 20; ++i)
        l_ordering.push_back(i);
    
    for(unsigned int i = 0; i < num_meioses; ++i)
        m_ordering.push_back(i);
    
    
    fprintf(stderr, "P(l-sampler) = %f\n", options.lsampler_prob);
    
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
                for(unsigned int k = l_ordering[j]; k < map->num_markers(); k += 20) {
                    lsamplers[k]->step(dg, k);
                }
            }
            /*
            // run every other locus independently
            for(unsigned int j = 0; j < 2; ++j) {
                #pragma omp parallel for
                for(unsigned int k = j; k < map->num_markers(); k += 2) {
                    lsamplers[k]->step(dg, k);
                }
            }
            */
        }
        else {
            random_shuffle(m_ordering.begin(), m_ordering.end());
            
            msampler.reset(dg, m_ordering[0]);
            for(unsigned int j = 0; j < m_ordering.size(); ++j) {
                msampler.step(dg, m_ordering[j]);
            }
        }
        
        p.increment();
        
        if(dg.get_likelihood() == LOG_ILLEGAL) {
            fprintf(stderr, "error: descent graph illegal...\n");
            abort();
        }
        
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

