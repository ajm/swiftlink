#ifndef LKG_MARKOVCHAIN_H_
#define LKG_MARKOVCHAIN_H_

using namespace std;

#include <cmath>
#include <cstdlib>


class SimwalkDescentGraph;
class Pedigree;
class GeneticMap;


class MarkovChain {

	Pedigree* ped;
	GeneticMap* map;
    double burnin_proportion;

    bool accept_metropolis(double new_prob, double old_prob);

 public:
	MarkovChain(Pedigree* p, GeneticMap* m) 
	    : ped(p), map(m), burnin_proportion(0.1) {}
    
	~MarkovChain() {}

    SimwalkDescentGraph* run(SimwalkDescentGraph* seed, unsigned iterations);
};

#endif

