#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

#include <time.h>

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
#include "omp_facade.h"

#ifdef USE_CUDA
  #include "gpu_lodscores.h"
#endif

using namespace std;


//#define CODA_OUTPUT 1

void MarkovChain::_init() {
    // heat up the map
    map.set_temperature(temperature);

    // create samplers
    for(int i = 0; i < min(get_max_threads(), int(map.num_markers())); ++i) {
        LocusSampler* tmp = new LocusSampler(ped, &map, psg, i, options.sex_linked);
        lsamplers.push_back(tmp);
    }

    lod = new LODscores(&map);

    // lod scorers
    for(int i = 0; i < min(get_max_threads(), int((map.num_markers() - 1) * map.get_lodscore_count())); ++i) {
        Peeler* tmp = new Peeler(ped, &map, psg, lod, options.sex_linked);
        peelers.push_back(tmp);
    }
    
    double trait_prob = peelers[0]->calc_trait_prob();
    
    printf("P(T) = %.5f\n", trait_prob / log(10));
    
    // lod score result objects
#ifdef USE_CUDA
    if(options.use_gpu) {
        gpulod = new GPULodscores(ped, &map, psg, options, trait_prob);
    }
#endif
    //lod = new LODscores(&map);
    lod->set_trait_prob(trait_prob);
    
    // ordering for locus samplers
    for(int i = 0; i < int(map.num_markers()); ++i) {
        l_ordering.push_back(i);
    }
    
    // ordering for meiosis samplers
    unsigned num_meioses = 2 * (ped->num_members() - ped->num_founders());
    for(unsigned int i = 0; i < num_meioses; ++i) {
        unsigned person_id = ped->num_founders() + (i / 2);
        enum parentage p = static_cast<enum parentage>(i % 2);
        
        Person* tmp = ped->get_by_index(person_id);
        
        if(not tmp->safe_to_ignore_meiosis(p)) {
            m_ordering.push_back(i);
        }
    }

    // create coda file if necessary
    if(options.coda_logging) {
        //string fname = options.coda_prefix + ".ped" + ped->get_id() + ".run" + to_string(seq_num); // does not work on gcc 4.8.X
        char buf[4];
        sprintf(buf, "%d", seq_num);
        string fname = options.coda_prefix + ".ped" + ped->get_id() + ".run" + string(buf);
        coda_filehandle = fopen(fname.c_str(), "w");
        fprintf(coda_filehandle, "iteration likelihood\n");
        printf("opened trace file (%s)\n", fname.c_str());
    }
}

void MarkovChain::_kill() {
    for(int i = 0; i < int(lsamplers.size()); ++i) {
        delete lsamplers[i];
    }
    for(int i = 0; i < int(peelers.size()); ++i) {
        delete peelers[i];
    }

    if(options.coda_logging) {
        fclose(coda_filehandle);
    }
}

void MarkovChain::step(DescentGraph& dg, int start_iteration, int step_size) {
    int thread_num = 0;

    for(int i = start_iteration; i < start_iteration + step_size; ++i) {
        if(get_random() < options.lsampler_prob) {
            
            random_shuffle(l_ordering.begin(), l_ordering.end());
            vector<int> thread_assignments(get_max_threads(), -1);
            vector<int> tmp(l_ordering);

            #pragma omp parallel
            {
                thread_num = get_thread_num(); 

                while(1) {
                    #pragma omp critical 
                    {
                        int locus = -1;
                        if(not tmp.empty()) {
                            for(int j = int(tmp.size()-1); j >= 0; --j) {
                                if(noninterferring(thread_assignments, tmp[j], thread_num)) {
                                    locus = tmp[j];
                                    tmp.erase(tmp.begin() + j);
                                    break;
                                }
                            }
                        }

                        thread_assignments[thread_num] = locus;
                    }
                    
                    if(thread_assignments[thread_num] == -1)
                        break;
                    
                    lsamplers[thread_num]->set_locus_minimal(thread_assignments[thread_num]);
                    lsamplers[thread_num]->step(dg, thread_assignments[thread_num]);
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
        
/*
#ifdef CODA_OUTPUT
        if((i % options.scoring_period) == 0) {
            double current_likelihood = dg.get_likelihood();
            if(current_likelihood == LOG_ILLEGAL) {
                fprintf(stderr, "error: descent graph illegal...\n");
                abort();
            }
            
            fprintf(stderr, "%d\t%f\n", i+1, current_likelihood);
        }

        continue;
#endif
*/
        if(i < options.burnin)
            continue;

        if(temperature != 1.0)
            continue;
        
        // only score the coldest chain
        if((i % options.scoring_period) == 0) {
         
#ifdef USE_CUDA
            if(not options.use_gpu) {
#endif
                #pragma omp parallel 
                {
                    thread_num = get_thread_num();
                    #pragma omp for
                    for(int j = 0; j < int(map.num_markers() - 1); ++j) {
                        peelers[thread_num]->set_locus(j);
                        peelers[thread_num]->process(&dg);
                    }
                }
#ifdef USE_CUDA
            }
            else {
                gpulod->calculate(dg);
            }
#endif
        }
        
    }
    
    
#ifdef USE_CUDA
    if(options.use_gpu) {
        gpulod->get_results(lod);
    }
#endif
}

void MarkovChain::run_scalable_lsampler(DescentGraph& dg, vector<int>& lgroups, int num_lgroups) {
    int thread_num = 0;

    random_shuffle(lgroups.begin(), lgroups.end());

    for(int j = 0; j < num_lgroups; ++j) {
        #pragma omp parallel num_threads(lsamplers.size()) private(thread_num)
        {   
            thread_num = get_thread_num();

            #pragma omp for
            for(int k = lgroups[j]; k < int(map.num_markers()); k += num_lgroups) {
                lsamplers[thread_num]->set_locus_minimal(l_ordering[k]);
                lsamplers[thread_num]->step(dg, l_ordering[k]);
            }
        }
    }
}

void MarkovChain::run_old_lsampler(DescentGraph& dg) {
    int thread_num = 0;
    random_shuffle(l_ordering.begin(), l_ordering.end());
    vector<int> thread_assignments(lsamplers.size(), -1);
    vector<int> tmp(l_ordering);

    #pragma omp parallel num_threads(lsamplers.size()) private(thread_num)
    {
        thread_num = get_thread_num();

        while(1) {
            #pragma omp critical 
            {
                int locus = -1;
                if(not tmp.empty()) {
                    for(int j = int(tmp.size()-1); j >= 0; --j) {
                        if(noninterferring(thread_assignments, tmp[j], thread_num)) {
                            locus = tmp[j];
                            tmp.erase(tmp.begin() + j);
                            break;
                        }
                    }
                }

                thread_assignments[thread_num] = locus;
            }
                    
            if(thread_assignments[thread_num] == -1)
                break;
                    
            lsamplers[thread_num]->set_locus_minimal(thread_assignments[thread_num]);
            lsamplers[thread_num]->step(dg, thread_assignments[thread_num]);

            #pragma omp critical
            {
                thread_assignments[thread_num] = -1;
            }
        }
    }
}

int MarkovChain::optimal_num_lgroups(DescentGraph& dg) {
    int best_num_lgroups = -1;
    double best_time, start_time, run_time;
    int repetitions = 100;

    if(get_max_threads() == 1) {
        return -1;
    }

    //best_time = numeric_limits<double>::infinity();

    // test the old version first
    start_time = get_wtime();
    for(int k = 0; k < repetitions; ++k) {
        run_old_lsampler(dg);
    }
    best_time = (get_wtime() - start_time) / repetitions;
    
//    fprintf(stderr, "%d lgroups, %f seconds\n", -1, best_time);

    // test scalable version with different step sizes
    for(int i = 3; i < 11; ++i) {
        vector<int> tmp(i);
        for(int j = 0; j < i; ++j) {
            tmp[j] = j;
        }

        start_time = get_wtime();
        for(int k = 0; k < repetitions; ++k) {
            run_scalable_lsampler(dg, tmp, i);
        }
        run_time = (get_wtime() - start_time) / repetitions;

        if(run_time < best_time) {
            best_num_lgroups = i;
            best_time = run_time;
        }

//        fprintf(stderr, "%d lgroups, %f seconds\n", i, run_time);
    }

    return best_num_lgroups;
}

// old version
LODscores* MarkovChain::run(DescentGraph& dg) {
    int thread_num = 0;
    int num_lgroups = -1;

    Progress p("MCMC: ", options.iterations + options.burnin);
    
    if(dg.get_likelihood() == LOG_ILLEGAL) {
        fprintf(stderr, "error: descent graph illegal pre-markov chain...\n");
        abort();
    }

//    fprintf(stderr, "\n");
    num_lgroups = optimal_num_lgroups(dg);
//    fprintf(stderr, "\noptimal number of lgroups = %d\n", num_lgroups);
//    exit(1);

    vector<int> lgroups; //(num_lgroups, 0);
    for(int i = 0; i < num_lgroups; ++i) {
        lgroups.push_back(i);
    }

    for(int i = 0; i < (options.iterations + options.burnin); ++i) {
        if(get_random() < options.lsampler_prob) {

            if(num_lgroups == -1) {
                run_old_lsampler(dg);
            }
            else {
                run_scalable_lsampler(dg, lgroups, num_lgroups);
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
        
        if(i < options.burnin) {
            continue;
        }
        
        if((i % options.scoring_period) == 0) {
            if(options.coda_logging) {
                double current_likelihood = dg.get_likelihood();
                if(current_likelihood == LOG_ILLEGAL) {
                    fprintf(stderr, "error: descent graph illegal...\n");
                    abort();
                }

                fprintf(coda_filehandle, "%d\t%f\n", i+1, current_likelihood);
            }
         
#ifdef USE_CUDA
            if(not options.use_gpu) {
#endif
                #pragma omp parallel private(thread_num)
                {
                    thread_num = get_thread_num();
                    #pragma omp for
                    for(int j = 0; j < int(map.num_markers() - 1); ++j) {
                        peelers[thread_num]->set_locus(j);
                        peelers[thread_num]->process(&dg);
                    }
                }
#ifdef USE_CUDA
            }
            else {
                gpulod->calculate(dg);
            }
#endif
        }
        
    }
    
    
    p.finish();
    
#ifdef USE_CUDA
    if(options.use_gpu) {
        gpulod->get_results(lod);
    }
#endif
    
    return lod;
}

bool MarkovChain::noninterferring(vector<int>& x, int val, int thread_num) {
    for(int i = 0; i < int(x.size()); ++i) {
        if((i != thread_num) and (x[i] != -1)) {
            int diff = val - x[i];
            if((diff == 1) or (diff == -1)) {
                return false;
            }
        }
    }
    return true;
}

