#ifndef LKG_GENETICALGORITHM_H_
#define LKG_GENETICALGORITHM_H_

#include <vector>

const double TRUNC_PROP = 0.1;

#include "descentgraph.h"

//class DescentGraph;
class GeneticMap;
class Pedigree;

// TODO: add a means to choose what selection method
// and crossover methods are called

// descent graph should have some kind of common interface superclass
// these methods should be in terms of that superclass
// or alternately could this be a template?

class GeneticAlgorithm {
    
    Pedigree* ped;
    GeneticMap* map;
    //unsigned popsize;
    unsigned truncation;
    
    void crossover_single(DescentGraph* p1, DescentGraph* p2, 
                          DescentGraph* c1, DescentGraph* c2);
    void crossover_double(DescentGraph* p1, DescentGraph* p2, 
                          DescentGraph* c1, DescentGraph* c2);
    void crossover_uniform(DescentGraph* p1, DescentGraph* p2, 
                           DescentGraph* c1, DescentGraph* c2);
    DescentGraph* truncation_selection(vector<DescentGraph>& pop);
	DescentGraph* fitness_proportionate_selection(vector<DescentGraph>& pop, 
														vector<double>& cdf);
    DescentGraph* deterministic_tournament_selection(vector<DescentGraph>& pop, 
														unsigned t);
	DescentGraph* probabilistic_tournament_selection(vector<DescentGraph>& pop, 
														unsigned t, double p);

    DescentGraph* mutate(DescentGraph* c, double prob);
    DescentGraph* mutate(DescentGraph* c, int freq);
	void make_cdf(vector<DescentGraph>& pop, vector<double>& cdf);
    unsigned get_random_nonfounder();
    unsigned get_random_locus();
    unsigned get_random(unsigned i);
    enum parentage get_random_parent();
    
 public :
    GeneticAlgorithm(Pedigree* p, GeneticMap* m) : ped(p), map(m) {}
    ~GeneticAlgorithm() {}
    
    //DescentGraph* optimise(unsigned popsize, unsigned iterations, 
    //                       unsigned elitism = 0);
    DescentGraph* optimise(unsigned population, unsigned iterations, unsigned elitism, 
		                  int crossover, int selection, unsigned selection_size, double selection_prob);
};

#endif
