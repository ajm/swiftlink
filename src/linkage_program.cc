using namespace std;

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "peel_sequence_generator.h"
#include "simwalk_descent_graph.h"
#include "simulated_annealing.h"
#include "linkage_program.h"
#include "markov_chain.h"
#include "genetic_map.h"
#include "peel_matrix.h"
#include "rfunction.h"
#include "pedigree.h"
#include "progress.h"
#include "peeling.h"


bool LinkageProgram::run() {
    bool ret = true;

    // XXX need to know how to do this properly, 
    // look up better random numbers for simulations etc
    srandom(time(NULL));

    for(unsigned int i = 0; i < pedigrees.size(); ++i) {
        ret &= run_pedigree(pedigrees[i]);
    }

    return ret;
}

bool LinkageProgram::run_pedigree(Pedigree& p) {
    unsigned int iterations = 10000; //800 * p.num_members() * p.num_markers() * 10 * 2;

    fprintf(stderr, "processing pedigree %s\n", p.get_id().c_str());

    // RUN SIMULATED ANNEALING
    SimulatedAnnealing sa(&p, &map);
    SimwalkDescentGraph* sdg1 = sa.optimise(iterations);

    printf("sa final prob = %f\n", sdg1->get_prob());
    
    // RUN MARKOV CHAIN
    MarkovChain mc(&p, &map);
    SimwalkDescentGraph* sdg2 = mc.run(sdg1, iterations);
    
    printf("mcmc final prob = %f\n", sdg2->get_prob());
    
    
    delete sdg1;

    
    // XXX put all peeling sequence calculation + evaluation code 
    // in own module, this is not the final form as it need to peel
    // over many descent graphs + all loci not just this one, but this is just
    // placeholder to get started...


    // PEELING + CALCULATION OF LOD SCORES
    PeelSequenceGenerator psg(p);
    psg.build_peel_order();
    psg.print();
        
    // setup r-functions
    vector<PeelOperation>& ops = psg.get_peel_order();
    
    
    vector<Rfunction> rfunctions;
    for(vector<PeelOperation>::size_type j = 0; j < ops.size(); ++j) {
        Rfunction rf(ops[j], &p, 4);
        rfunctions.push_back(rf);
    }

    // perform the peel for every locus
    Progress prog("Peeling:", rfunctions.size());
    PeelMatrix* last = NULL;
    for(vector<Rfunction>::size_type j = 0; j < rfunctions.size(); ++j) {
        Rfunction& rf = rfunctions[j];
            
        if(not rf.evaluate(last)) {
            prog.error("bad R function");
            goto peel_error;
        }
        
        prog.increment();
        
        last = rf.get_matrix();
    }

    prog.finish();


    delete sdg2;
	return true;


peel_error :

    delete sdg2;
    return false;
}

