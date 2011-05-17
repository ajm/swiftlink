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
    double r = log(random() / double(RAND_MAX));

    return r < (new_prob - old_prob);
}

Peeler* MarkovChain::run(DescentGraph* seed, unsigned iterations) {
	SimwalkSampler ss(ped, map);
	DescentGraph current(*seed);
	DescentGraphDiff dgd;
	double prob;
    unsigned burnin_steps = unsigned(iterations * burnin_proportion); // XXX 
    

    Progress p("Markov Chain:", iterations);
    p.start();

	for(unsigned i = 0; i < iterations; ++i) {

        ss.step(dgd);
        
        p.increment();
        
        if((i % 1000) == 0) {
            peel.process(&current); // XXX use reference?
        }
        
        if(not current.evaluate_diff(dgd, &prob)) {
            continue;
        }
        
        // separate for now, just in case I want to add any stat gathering...
		if(i < burnin_steps) {
		    current.apply_diff(dgd);
		    /*
		    double tmp_prob;
		    current.likelihood(&tmp_prob);
            printf("%d: complete = %f, diff = %f\n", i, tmp_prob, prob);
            */
		    continue;
        }
        
        if(accept_metropolis(prob, current.get_prob())) {
            current.apply_diff(dgd);
        }
        /*
        if((i % 1000) == 0) {
            peel.process(&current); // XXX use reference?
        }
        */
	}

    p.finish();
	
	return &peel;
}

