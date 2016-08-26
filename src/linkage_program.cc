#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>

#include "descent_graph.h"
#include "linkage_program.h"
#include "linkage_writer.h"
#include "disease_model.h"
#include "markov_chain.h"
#include "genetic_map.h"
#include "pedigree.h"
#include "peeler.h"
#include "random.h"
#include "sequential_imputation.h"
//#include "gpu_markov_chain.h"
#include "lod_score.h"

#include "mc3.h"

using namespace std;


bool LinkageProgram::run() {
    vector<LODscores*> all_scores;
    LODscores* tmp;
    bool ret = true;
    LinkageWriter lw(&map, outfile, options.verbose);


    fprintf(stderr, "\nLinkage parameters:\n"
                    "\tpenetrance = %.2f:%.2f:%.2f\n"
                    "\ttrait freq = %.2e\n"
                    "\tsex-linked = %s\n"
                    "\tburnin iterations = %d\n"
                    "\tsampling iterations = %d\n"
                    "\tsampling period = %d\n"
                    "\tlocus sampler prob = %.3f\n"
                    "\tnumber of runs = %d\n\n",
                    dm.get_penetrance(TRAIT_HOMO_U),
                    dm.get_penetrance(TRAIT_HETERO),
                    dm.get_penetrance(TRAIT_HOMO_A),
                    dm.get_freq(),
                    options.sex_linked ? "true" : "false",
                    options.burnin,
                    options.iterations,
                    options.scoring_period,
                    options.lsampler_prob,
                    options.mcmc_runs);


    init_random();
    if(options.random_filename == "") {
        seed_random_implicit();
    }
    else {
        seed_random_explicit(options.random_filename);
    }
    
    for(unsigned int i = 0; i < pedigrees.size(); ++i) {
        
        if(options.verbose) {
            fprintf(stderr, "%s\n", pedigrees[i].debug_string().c_str());
        }
        
        // it cannot actually be NULL, the program will call
        // abort() at the slightest hint of a problem
        if((tmp = run_pedigree_average(pedigrees[i], options.mcmc_runs)) == NULL) {
            fprintf(stderr, "error: pedigree '%s' failed\n", pedigrees[i].get_id().c_str());
            
            ret = false;
            goto die;
        }
        
        all_scores.push_back(tmp);
    }
    
    //LinkageWriter lw(&map, outfile, options.verbose);
    if(not lw.write(all_scores)) {
        fprintf(stderr, "error: could not write output file '%s'\n", outfile.c_str());
        
        ret = false;
        goto die;
    }
    
die:
    for(unsigned int i = 0; i < all_scores.size(); ++i) {
        delete all_scores[i];
    }
    
    return ret;
}

LODscores* LinkageProgram::run_pedigree_average(Pedigree& p, int repeats) {
    LODscores *ret, *tmp;

    ret = run_pedigree(p, 0);

    for(int i = 1; i < repeats; ++i) {
        tmp = run_pedigree(p, i);
        ret->merge_results(tmp);
        delete tmp;
    }

    return ret;
}

LODscores* LinkageProgram::run_pedigree(Pedigree& p, int sequence_number) {
    
    if(options.verbose) {
        fprintf(stderr, "processing pedigree %s\n", p.get_id().c_str());
    }

    if(options.affected_only) {
        for(unsigned int i = 0; i < p.num_members(); ++i) {
            Person* q = p.get_by_index(i);

            if(not q->isaffected()) {
                q->make_unknown_affection();
            }
        }
    }


    DescentGraph dg(&p, &map, dm.is_sexlinked());
    dg.random_descentgraph(); // just in case the user selects zero sequential imputation iterations
    
    if(dg.get_likelihood() == LOG_ZERO) {
        fprintf(stderr, "error, bad descent graph %s:%d\n", __FILE__, __LINE__);
        abort();
    }

    
    PeelSequenceGenerator psg(&p, &map, dm.is_sexlinked(), options.verbose);
    psg.build_peel_sequence(options.peelopt_iterations);

    if(options.verbose) {
        fprintf(stderr, "\n\n%s\n\n", psg.debug_string().c_str());
    }


    SequentialImputation si(&p, &map, &psg, dm.is_sexlinked());
    //si.run(dg, options.si_iterations);
    si.parallel_run(dg, options.si_iterations);

    if(dg.get_likelihood() == LOG_ZERO) {
        fprintf(stderr, "error, bad descent graph %s:%d\n", __FILE__, __LINE__);
        abort();
    }


    /*
    if(not options.use_gpu) {
        MarkovChain chain(&p, &map, &psg, options);
        return chain.run(dg);
    }
    else {
        GPUMarkovChain chain(&p, &map, &psg, options);
        return chain.run(dg);
    }
    
    fprintf(stderr, "error: nothing was run (%s:%d)\n", __FILE__, __LINE__);
    abort();
    */
    
    options.sex_linked = dm.is_sexlinked();

    MarkovChain chain(&p, &map, &psg, options, sequence_number);
    return chain.run(dg);

    //Mc3 chain(&p, &map, &psg, options);
    //return chain.run();
}

