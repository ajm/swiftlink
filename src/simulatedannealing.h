#ifndef LKG_SIMULATEDANNEALING_H_
#define LKG_SIMULATEDANNEALING_H_

using namespace std;

#include <cmath>

#include "SA_descentgraph.h"
#include "pedigree.h"
#include "geneticmap.h"


const double START_TEMPERATURE = 100.0;
const double TEMPERATURE_CHANGE_FACTOR = 0.99;
const int TEMPERATURE_CHANGES = 800;

class SimulatedAnnealing {

	Pedigree* ped;
	GeneticMap* map;
	double temperature;

 public:
	SimulatedAnnealing(Pedigree* p, GeneticMap* m) 
	    : ped(p), map(m), temperature(START_TEMPERATURE) {}
	~SimulatedAnnealing() {}

	bool accept_metropolis(double new_prob, double old_prob) {
		return log(rand() / double(RAND_MAX)) < (new_prob - old_prob);
	}
	
    bool accept_annealing(double new_prob, double old_prob, double temp) {
        if(old_prob < new_prob)
            return true;
        
        double r = log(rand() / double(RAND_MAX)) * temp;
        double p = new_prob - old_prob;
        
		return r < p;
	}

	SA_DescentGraph* optimise(unsigned iterations) {
		SA_DescentGraph* best;
		SA_DescentGraph* current;
		SA_DescentGraph* temp;
		double prob;

		current = new SA_DescentGraph(ped, map);
		current->random();
		current->haplotype_likelihood();

		temp = new SA_DescentGraph(ped, map);
		best = new SA_DescentGraph(ped, map);
		*best = *current;

		for(unsigned i = 0; i < iterations; ++i) {
			*temp = *current;
			
			if((i % (iterations / TEMPERATURE_CHANGES)) == 0) {
			    temperature *= TEMPERATURE_CHANGE_FACTOR;
            }

			temp->step();
			temp->haplotype_likelihood(&prob); 
			// TODO temporary, put this in SA_DescentGraph only on affected loci

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
			
			if((i % 10000) == 0) {
			    fprintf(stdout, "%d %f %f (%f)\n", \
			                i, \
			                current->get_prob() / log(10), \
			                best->get_prob() / log(10), \
			                temperature);
			}
		}
		
		delete temp;
		delete current;

		return best;
	}
};

#endif

