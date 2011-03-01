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
#include "peeler.h"


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

    printf("processing pedigree %s\n", p.get_id().c_str());

    // RUN SIMULATED ANNEALING
    SimulatedAnnealing sa(&p, &map);
    SimwalkDescentGraph* sdg1 = sa.optimise(iterations);

    printf("sa final prob = %f\n", sdg1->get_prob());
    
    // RUN MARKOV CHAIN
    MarkovChain mc(&p, &map);
    SimwalkDescentGraph* sdg2 = mc.run(sdg1, iterations);
    
    printf("mcmc final prob = %f\n", sdg2->get_prob());
    
    delete sdg1;

    Peeler peeler(&p, &map);
    if(not peeler.peel(sdg2)) {
        delete sdg2;    
        return false;
    }

    delete sdg2;
	return true;
}

