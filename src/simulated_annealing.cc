using namespace std;

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "pedigree.h"
#include "genetic_map.h"
#include "simwalk_descent_graph.h"
#include "simulated_annealing.h"


bool SimulatedAnnealing::accept_metropolis(double new_prob, double old_prob) {
	return log(random() / double(RAND_MAX)) < (new_prob - old_prob);
}
	
bool SimulatedAnnealing::accept_annealing(double new_prob, double old_prob, double temp) {
    if(old_prob < new_prob)
        return true;
    
    double r = log(random() / double(RAND_MAX)) * temp;
    double p = new_prob - old_prob;
    
    return r < p;
}

SimwalkDescentGraph* SimulatedAnnealing::optimise(unsigned iterations) {
    SimwalkDescentGraph* best;
	SimwalkDescentGraph* current;
	SimwalkDescentGraph* temp;
	double prob;

	current = new SimwalkDescentGraph(ped, map);
	current->random();
	current->haplotype_likelihood();

	temp = new SimwalkDescentGraph(ped, map);
	best = new SimwalkDescentGraph(ped, map);

	*best = *current;

	for(unsigned i = 0; i < iterations; ++i) {
		*temp = *current;
			
		if((i % (iterations / TEMPERATURE_CHANGES)) == 0) {
		    temperature *= TEMPERATURE_CHANGE_FACTOR;
        }

		temp->step();
		//temp->haplotype_likelihood(&prob);
        temp->likelihood(&prob); 
		
// XXX of course, this all involves massive amounts of copying
// need everything to work on a per-loci basis, ie: likelihood calculations
// and all this...
		if(temp->illegal()) {
			continue;
		}
			
		if(accept_annealing(temp->get_prob(), current->get_prob(), temperature)) {
			*current = *temp;
/*
			fprintf(stdout, "%d %f %f (%f) ACCEPT\n", \
		                i, \
		                current->get_prob() / log(10), \
		                best->get_prob() / log(10), \
		                temperature);
*/
        }

		if(best->get_prob() < current->get_prob()) {
		    *best = *current;
		}
/*
		if((i % 10000) == 0) {
		    fprintf(stdout, "%d %f %f (%f)\n", \
		                i, \
		                current->get_prob() / log(10), \
		                best->get_prob() / log(10), \
		                temperature);
    	}
*/
	}
		
	delete temp;
	delete current;

	return best;
}

