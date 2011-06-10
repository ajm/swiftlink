#ifndef LKG_MARKOVCHAIN_H_
#define LKG_MARKOVCHAIN_H_

using namespace std;

#include <cmath>
#include <cstdlib>

#include "peeler.h"

class DescentGraph;
class Pedigree;
class GeneticMap;

const unsigned int SAMPLE_PERIOD = 1000;


class MarkovChain {

	Pedigree& ped;
	GeneticMap& map;
	Peeler peel;
    double burnin_proportion;
    unsigned illegal;
    unsigned accepted;
    unsigned rejected;

    bool accept_metropolis(double new_prob, double old_prob);

 public:
	MarkovChain(Pedigree& p, GeneticMap& m) 
	    : ped(p), map(m), peel(p, m), burnin_proportion(0.2), 
	      illegal(0), accepted(0), rejected(0) {}
    
	~MarkovChain() {}

    Peeler* run(DescentGraph* seed, unsigned iterations);
};

#endif

