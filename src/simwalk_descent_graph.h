#ifndef LKG_SA_DESCENTGRAPH_H_
#define LKG_SA_DESCENTGRAPH_H_

using namespace std;

#include <cstdlib>
#include <vector>

#include "genotype.h"
#include "pedigree.h"
#include "person.h"
#include "descent_graph.h"


class SimwalkDescentGraph : public DescentGraph {

    unsigned _geometric();
    int select_transition_rule();
    inline unsigned get_random(unsigned i);
    unsigned get_random_person();
    unsigned get_random_nonleaf();
    unsigned get_random_nonfounder();
    unsigned get_random_locus();
    enum parentage get_random_parent();void transition_t0(unsigned locus);
    void transition_t0(unsigned id, unsigned locus, enum parentage p);
    void transition_t1(unsigned locus);
    void transition_t1(Person* p, unsigned locus);
    void transition_t2a(unsigned locus);
    void transition_t2b(unsigned locus);
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

