#ifndef LKG_SIMULATEDANNEALING_H_
#define LKG_SIMULATEDANNEALING_H_

using namespace std;

#include <cmath>
#include <cstdlib>


class DescentGraph;
class Pedigree;
class GeneticMap;

const double START_TEMPERATURE = 100.0;
const double TEMPERATURE_CHANGE_FACTOR = 0.99;
const unsigned int TEMPERATURE_CHANGES = 800;


class SimulatedAnnealing {

	Pedigree* ped;
	GeneticMap* map;
	double temperature;

    bool accept_metropolis(double new_prob, double old_prob);
    bool accept_annealing(double new_prob, double old_prob, double temp);

 public:
	SimulatedAnnealing(Pedigree* p, GeneticMap* m) : 
	    ped(p), 
	    map(m), 
	    temperature(START_TEMPERATURE) {}
    
    SimulatedAnnealing(const SimulatedAnnealing& rhs) :
        ped(rhs.ped),
        map(rhs.map),
        temperature(rhs.temperature) {}
    
	~SimulatedAnnealing() {}
	
	SimulatedAnnealing& operator=(const SimulatedAnnealing& rhs) {
	    
	    if(&rhs != this) {
	        ped = rhs.ped;
	        map = rhs.map;
	        temperature = rhs.temperature;
	    }
	    
	    return *this;
	}

    DescentGraph* optimise(unsigned iterations);
};

#endif

