#ifndef LKG_SIMWALKSAMPLER_H_
#define LKG_SIMWALKSAMPLER_H_

using namespace std;

#include <cstdlib>
#include <vector>

#include "genotype.h"
#include "pedigree.h"
#include "person.h"
#include "genetic_map.h"
#include "descent_graph_diff.h"


class SimwalkSampler {

    Pedigree* ped;
    GeneticMap* map;

    unsigned _geometric();
    unsigned select_next_locus(unsigned locus);
    int select_transition_rule(unsigned person);
    inline unsigned get_random(unsigned i);
    unsigned get_random_person();
    unsigned get_random_nonfounder();
    unsigned get_random_locus();
    enum parentage get_random_parent();
    DescentGraphDiff transition_t0(unsigned person, unsigned locus);
    void transition_t0(unsigned id, unsigned locus, enum parentage p);
/*
    void transition_t1(unsigned person, unsigned locus);
    void transition_t1(Person* p, unsigned locus);
    void transition_t2a(unsigned person, unsigned locus);
    void transition_t2b(unsigned person, unsigned locus);
    void transition_t2(unsigned id, unsigned locus, bool same_gender);
*/

 public:
	SimwalkSampler(Pedigree* ped, GeneticMap* map) 
		: ped(ped), map(map) {}
    
	~SimwalkSampler() {}
	
    DescentGraphDiff step();
};

#endif

