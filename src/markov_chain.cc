using namespace std;

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "pedigree.h"
#include "genetic_map.h"
#include "simwalk_descent_graph.h"
#include "markov_chain.h"
#include "progress.h"


bool MarkovChain::accept_metropolis(double new_prob, double old_prob) {
	return log(random() / double(RAND_MAX)) < (new_prob - old_prob);
}

SimwalkDescentGraph* MarkovChain::run(SimwalkDescentGraph* seed, unsigned iterations) {
    SimwalkDescentGraph* best;
	SimwalkDescentGraph* current;
	SimwalkDescentGraph* temp;
	double prob;
    unsigned burnin_steps = iterations * burnin_proportion;
    
    
    current = new SimwalkDescentGraph(ped, map);
    *current = *seed;

	temp = new SimwalkDescentGraph(ped, map);
	best = new SimwalkDescentGraph(ped, map);

	*best = *current;

    Progress p("Markov Chain:", iterations);
    p.start();

	for(unsigned i = 0; i < iterations; ++i) {
		*temp = *current;

		temp->step();
        temp->likelihood(&prob); 
		
        p.increment();

        // XXX of course, this all involves massive amounts of copying
        // need everything to work on a per-loci basis, ie: likelihood calculations
        // and all this...
		if(temp->illegal()) {
			continue;
		}
	
        // separate for now, just in case I want to add any stat gathering...
		if(i < burnin_steps) {
			*current = *temp;
        }
        else if(accept_metropolis(temp->get_prob(), current->get_prob())) {
            *current = *temp;
        }

		if(best->get_prob() < current->get_prob()) {
		    *best = *current;
		}
	}

    p.finish();
		
	delete temp;
	delete current;

	return best;
}

