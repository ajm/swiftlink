#ifndef LKG_SIMWALKDESCENTGRAPH_H_
#define LKG_SIMWALKDESCENTGRAPH_H_

using namespace std;

#include <cstdlib>
#include <vector>

#include "genotype.h"
#include "pedigree.h"
#include "person.h"
#include "descent_graph.h"


class SimwalkDescentGraph : public DescentGraph {

    unsigned _geometric();
    unsigned select_next_locus(unsigned locus);
    int select_transition_rule(unsigned person);
    inline unsigned get_random(unsigned i);
    unsigned get_random_person();
    unsigned get_random_nonleaf();
    unsigned get_random_nonfounder();
    unsigned get_random_locus();
    enum parentage get_random_parent();
    void transition_t0(unsigned person, unsigned locus);
    void transition_t0(unsigned id, unsigned locus, enum parentage p);
    void transition_t1(unsigned person, unsigned locus);
    void transition_t1(Person* p, unsigned locus);
    void transition_t2a(unsigned person, unsigned locus);
    void transition_t2b(unsigned person, unsigned locus);
    void transition_t2(unsigned id, unsigned locus, bool same_gender);

 public:
	SimwalkDescentGraph(Pedigree* ped, GeneticMap* map) 
		: DescentGraph(ped, map) {}
    
	SimwalkDescentGraph(const SimwalkDescentGraph& sdg) 
		: DescentGraph(sdg) {}

	~SimwalkDescentGraph() {}
	
	SimwalkDescentGraph& operator=(const SimwalkDescentGraph& rhs) {
	    //printf("SA_DescentGraph=\n");
		if(this != &rhs) {
			DescentGraph::operator=(rhs);
		}
        
		return *this;
	}
	
    void step();
};

#endif

