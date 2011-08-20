#ifndef LKG_SIMWALKSAMPLER_H_
#define LKG_SIMWALKSAMPLER_H_

using namespace std;

#include <cstdlib>
#include <vector>

#include "genotype.h"
#include "pedigree.h"
#include "person.h"
#include "descent_graph_diff.h"
#include "descent_graph_types.h"


class SimwalkSampler {

    Pedigree* ped;
    DescentGraph* dg;
    
    unsigned type0_count;
    unsigned type1_count;
    unsigned type2_count;

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
    void transition_t2(DescentGraphDiff& dgd, unsigned person, unsigned locus, bool same_gender);

 public:
	SimwalkSampler(Pedigree* ped, DescentGraph* dg) : 
	    ped(ped), 
		dg(dg), 
		type0_count(0), 
		type1_count(0), 
		type2_count(0) {}
		  
    SimwalkSampler(const SimwalkSampler& rhs) : 
        ped(rhs.ped), 
		dg(rhs.dg), 
		type0_count(rhs.type0_count), 
		type1_count(rhs.type1_count), 
		type2_count(rhs.type2_count) {}
    
	~SimwalkSampler() {}
	
	SimwalkSampler& operator=(const SimwalkSampler& rhs) {
	    
	    if(&rhs != this) {
	        ped = rhs.ped;
	        dg = rhs.dg;
	        type0_count = rhs.type0_count;
	        type1_count = rhs.type1_count;
	        type2_count = rhs.type2_count;
	    }
	    
	    return *this;
	}
	
    void step(DescentGraphDiff& dgd);
};

#endif

