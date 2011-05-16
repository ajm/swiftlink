using namespace std;

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "pedigree.h"
#include "genetic_map.h"
#include "descent_graph.h"
#include "descent_graph_diff.h"
#include "simulated_annealing.h"
#include "simwalk_sampler.h"
#include "progress.h"


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

DescentGraph* SimulatedAnnealing::optimise(unsigned iterations) {
	SimwalkSampler ss(ped, map);
	DescentGraph current(ped, map);
	DescentGraph* best;
	double prob;

	current.random_descentgraph();
	//current->haplotype_likelihood();
    current.likelihood();
    
	best = new DescentGraph(current);
    
    if(iterations < TEMPERATURE_CHANGES) {
        iterations = TEMPERATURE_CHANGES;
    }
    
    Progress p("Simulated Annealing:", iterations);
    p.start();

	for(unsigned i = 0; i < iterations; ++i) {
		
		if((i % (iterations / TEMPERATURE_CHANGES)) == 0) {
		    temperature *= TEMPERATURE_CHANGE_FACTOR;
        }

		DescentGraphDiff dgd = ss.step();
        
        p.increment();
        
        if(not current.evaluate_diff(dgd, &prob)) {
            continue;
        }
		
		if(accept_annealing(prob, current.get_prob(), temperature)) {
			current.apply_diff(dgd);
			/*
			double tmp_prob;
		    current.likelihood(&tmp_prob);
    		printf("%d: complete = %f, diff = %f\n", i, tmp_prob, prob);
		    */
        }
        
		if(best->get_prob() < current.get_prob()) {
		    *best = current;
		}
		
		// finish early if there are no recombinations
		// XXX could a descent graph exist that has a lower likelihood, but a 
		// few recombinations?
		if(current.num_recombinations() == 0) {
		    fprintf(stderr, "\nsimulated annealing, early termination\n");
		    break;
		}
		
		//fprintf(stderr, "%d %d\n", i, current.num_recombinations());
	}

    p.finish();

	return best;
}

