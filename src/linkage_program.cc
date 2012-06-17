using namespace std;

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
#include "gpu_markov_chain.h"


bool LinkageProgram::run() {
    vector<double*> lod_scores;
    double* tmp;
    
    /*
    if(options.verbose) {
        fprintf(stderr, "%s\n", dm.debug_string().c_str());
        fprintf(stderr, "%s\n", map.debug_string().c_str());
    }
    */
        
    init_random();
    seed_random_implicit();
    
    for(unsigned int i = 0; i < pedigrees.size(); ++i) {
        
        if(options.verbose) {
            fprintf(stderr, "%s\n", pedigrees[i].debug_string().c_str());
        }
        
        
        // it cannot actually be NULL, the program will call
        // abort() at the slightest hint of a problem
        if((tmp = run_pedigree(pedigrees[i])) == NULL) {
            fprintf(stderr, "error: pedigree '%s' failed\n", pedigrees[i].get_id().c_str());
            free_lodscores(lod_scores);
            
            return false;
        }
        
        lod_scores.push_back(tmp);
    }
    
    LinkageWriter lw(&map, outfile, options.verbose);
    if(not lw.write(lod_scores)) {
        fprintf(stderr, "error: could not write output file '%s'\n", outfile.c_str());
        
        return false;
    }
    
    free_lodscores(lod_scores);
    
    return true;
}

double* LinkageProgram::run_pedigree(Pedigree& p) {
    
    if(options.verbose) {
        fprintf(stderr, "processing pedigree %s\n", p.get_id().c_str());
    }
    
    DescentGraph dg(&p, &map);
    
    dg.random_descentgraph(); // just in case the user selects zero sequential imputation iterations
    
    if(dg.get_likelihood() == LOG_ZERO) {
        fprintf(stderr, "error, bad descent graph %s:%d\n", __FILE__, __LINE__);
        abort();
    }
    
    
    PeelSequenceGenerator psg(&p, &map, options.verbose);
    if((options.peelseq_filename == "") or not psg.read_from_file(options.peelseq_filename)) {
        psg.build_peel_order();
    }
    
    SequentialImputation si(&p, &map, &psg);
    //si.run(dg, options.si_iterations);
    si.parallel_run(dg, options.si_iterations);
    

    if(not options.use_gpu) {
        MarkovChain chain(&p, &map, &psg, options);
        return chain.run(dg);
    }
    else {
        //GPUMarkovChain chain(&p, &map, &psg, options);
        //return chain.run(dg);
    }
    
    fprintf(stderr, "error: nothing was run (%s:%d)\n", __FILE__, __LINE__);
    abort();
}

void LinkageProgram::free_lodscores(vector<double*>& x) {
    for(unsigned i = 0; i < x.size(); ++i) {
        delete[] x[i];
    }
}

