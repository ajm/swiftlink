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
    
    return (log(random() / double(RAND_MAX)) * temp) < (new_prob - old_prob);
}

DescentGraph* SimulatedAnnealing::optimise(unsigned iterations) {
	DescentGraph current(ped, map);
	DescentGraph* best;
	DescentGraphDiff dgd;
	SimwalkSampler ss(ped, &current);
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

		ss.step(dgd);
		
		//printf("\n\n");
		//dgd.print();
	    //exit(0);
        
        p.increment();
        
        if(not current.evaluate_diff(dgd, &prob)) {
            continue;
        }
		
		if(accept_annealing(prob, current.get_prob(), temperature)) {
			current.apply_diff(dgd);
        }
        
		if(best->get_prob() < current.get_prob()) {
		    *best = current;
		}
        
		// finish early if there are no recombinations
		if(best->num_recombinations() == 0) {
		    fprintf(stderr, "\nsimulated annealing, early termination\n");
		    break;
		}
	}

    p.finish();

	return best;
}

