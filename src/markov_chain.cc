using namespace std;

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "pedigree.h"
#include "genetic_map.h"
#include "descent_graph.h"
#include "descent_graph_diff.h"
#include "markov_chain.h"
#include "progress.h"
#include "peeler.h"
#include "simwalk_sampler.h"


bool MarkovChain::accept_metropolis(double new_prob, double old_prob) {
    return (log(random() / double(RAND_MAX))) < (new_prob - old_prob);
}

Peeler* MarkovChain::run(DescentGraph* seed, unsigned iterations) {
	DescentGraph current(*seed);
	DescentGraphDiff dgd;
	SimwalkSampler ss(ped, current);
	double prob;
    unsigned burnin_steps = unsigned(iterations * burnin_proportion); // XXX 
    

    Progress p("Markov Chain:", iterations);
    p.start();

	for(unsigned i = 0; i < iterations; ++i) {

        ss.step(dgd);
        
        p.increment();
        
        if((i >= burnin_steps) and ((i % SAMPLE_PERIOD) == 0)) {
            peel.process(current);
        }
        
        if(not current.evaluate_diff(dgd, &prob)) {
            continue;
        }
        
        // separate for now, just in case I want to add any stat gathering...
		if(i < burnin_steps) {
		    current.apply_diff(dgd);
		    continue;
        }
        
        if(accept_metropolis(prob, current.get_prob())) {
            current.apply_diff(dgd);
        }
	}

    p.finish();
	
	return &peel;
}

