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
/*
    fprintf(stderr, "\n");
    fprintf(stderr, "r = %e\n", r);
    fprintf(stderr, "dx = %e\n", new_prob - old_prob));
    fprintf(stderr, "%s\n", r < (new_prob - old_prob) ? "ACCEPT" : "REJECT");
*/
    return r < (new_prob - old_prob);

//	return log(random() / double(RAND_MAX)) < (new_prob - old_prob);
}

void MarkovChain::run(SimwalkDescentGraph* seed, unsigned iterations) {
//    SimwalkDescentGraph* best;
	SimwalkDescentGraph* current;
	SimwalkDescentGraph* temp;
	double prob;
    unsigned burnin_steps = unsigned(iterations * burnin_proportion);
    
    
    current = new SimwalkDescentGraph(ped, map);
    *current = *seed;

	temp = new SimwalkDescentGraph(ped, map);
//	best = new SimwalkDescentGraph(ped, map);

//	*best = *current;

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
		    peel.peel(current);
		}

        // XXX of course, this all involves massive amounts of copying
        // need everything to work on a per-loci basis, ie: likelihood calculations
        // and all this...
		if(temp->illegal()) {
//		    fprintf(stderr, "ILLEGAL\n");
			continue;
		}
		
//		fprintf(stderr, "LEGAL\n");

        //printf("%e\n", current->get_prob());

	
        // separate for now, just in case I want to add any stat gathering...
		if(i < burnin_steps) {
			*current = *temp;
        }
        else if(accept_metropolis(temp->get_prob(), current->get_prob())) {
            *current = *temp;
            fprintf(stderr, "\ncurrent prob = %f\n", current->get_prob());
            current->print();
        }
/*
		if(best->get_prob() < current->get_prob()) {
		    *best = *current;
		}
*/
	}

    p.finish();
		
	delete temp;
	delete current;

//	return best;
}

Peeler* MarkovChain::get_peeler() {
    return &peel;
}


