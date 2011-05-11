using namespace std;

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "pedigree.h"
#include "genetic_map.h"
#include "simwalk_descent_graph.h"
#include "markov_chain.h"
#include "progress.h"
#include "peeler.h"


bool MarkovChain::accept_metropolis(double new_prob, double old_prob) {
    double r = log(random() / double(RAND_MAX));

    return r < (new_prob - old_prob);
}

Peeler* MarkovChain::run(SimwalkDescentGraph* seed, unsigned iterations) {
	SimwalkDescentGraph* current;
	SimwalkDescentGraph* temp;
	double prob;
    unsigned burnin_steps = unsigned(iterations * burnin_proportion);
    
    
    current = new SimwalkDescentGraph(ped, map);
    *current = *seed;

	temp = new SimwalkDescentGraph(ped, map);


    Progress p("Markov Chain:", iterations);
    p.start();

	for(unsigned i = 0; i < iterations; ++i) {
		*temp = *current;

		temp->step();
        temp->likelihood(&prob); 
		
        p.increment();
        
        // perform peeling
		if((i >= burnin_steps) and ((i % 1000) == 0)) {
		    //printf("peel!\n");
		    peel.process(static_cast<DescentGraph*>(current));
		}

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
	}

    p.finish();
		
	delete temp;
	delete current;
	
	return &peel;
}

