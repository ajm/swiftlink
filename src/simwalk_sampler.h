#ifndef LKG_SIMWALKSAMPLER_H_
#define LKG_SIMWALKSAMPLER_H_

using namespace std;

#include <cstdlib>
#include <vector>

#include "genotype.h"
#include "pedigree.h"
#include "person.h"
//#include "genetic_map.h"
#include "descent_graph_diff.h"
#include "descent_graph_types.h"


class SimwalkSampler {

    Pedigree* ped;
    DescentGraph* dg;

    unsigned _geometric();
    unsigned select_next_locus(unsigned locus);
    int select_transition_rule(unsigned person);
    inline unsigned get_random(unsigned i);
    unsigned get_random_person();
    unsigned get_random_nonfounder();
    unsigned get_random_locus();
    enum parentage get_random_parent();
    
    void transition_t0(DescentGraphDiff& dgd, unsigned person, unsigned locus);
    void transition_t0(DescentGraphDiff& dgd, unsigned person, unsigned locus, enum parentage parent);
    void transition_t1(DescentGraphDiff& dgd, unsigned person, unsigned locus);
    void transition_t2a(DescentGraphDiff& dgd, unsigned person, unsigned locus);
    void transition_t2b(DescentGraphDiff& dgd, unsigned person, unsigned locus);
    void transition_t2(DescentGraphDiff& dgd, unsigned id, unsigned locus, bool same_gender);

 public:
	SimwalkSampler(Pedigree* ped, DescentGraph* dg) 
		: ped(ped), dg(dg) {}
    
	~SimwalkSampler() {}
	
    void step(DescentGraphDiff& dgd);
};

#endif

