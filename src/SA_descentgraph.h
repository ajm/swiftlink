#ifndef LKG_SA_DESCENTGRAPH_H_
#define LKG_SA_DESCENTGRAPH_H_

using namespace std;

#include <cstdlib>
#include <vector>

#include "genotype.h"
#include "pedigree.h"
#include "person.h"
#include "descentgraph.h"

class SA_DescentGraph : public DescentGraph {

 public:
	SA_DescentGraph(Pedigree* ped, GeneticMap* map) 
		: DescentGraph(ped, map) {}
	SA_DescentGraph(const SA_DescentGraph& sdg) 
		: DescentGraph(sdg) {}

	~SA_DescentGraph() {}
	
	SA_DescentGraph& operator=(const SA_DescentGraph& rhs) {
	    //printf("SA_DescentGraph=\n");
		if(this != &rhs) {
			DescentGraph::operator=(rhs);
		}
		return *this;
	}
	
	// go geom(0.5) transitions per step of the chain
	// select people + transition rules at random, but not necessarily uniform
	// 		perform t0 > t1 > t2
	// 		select untyped more than typed
	// person + locus need in all rules, add as parameters
	
	// t0 - non-founders (founders are a waste of time)
	// t1 - non-leaves
	// t2 - non-leaves (people with spouses)
	
	void step() {
		vector<unsigned> loci;
		unsigned steps = _geometric();
		
		//printf("num locus = %d\n", steps);
		
		//unsigned l = get_random_locus();
		for(unsigned i = 0; i < steps; ++i) {
			unsigned l = get_random_locus();

			switch(select_transition_rule()) {
				case 0 :
					transition_t0(l);
					break;
				case 1 :
					transition_t1(l);
					break;
				case 2 :
					transition_t2a(l);
					break;
				case 3 :
					transition_t2b(l);
					break;
				default:
					break;
			}

			loci.push_back(l);
/*
			if(l < ped->num_markers())
				l++;
			else 
				l = get_random_locus();
*/
		}
		
		// TODO recalculate likelihood for all loci
		//for(unsigned i = 0; i < loci.size(); ++i) {}
	}

 private:
	unsigned _geometric() {
		unsigned c = 1;

		while(1) {
			if((random() / double(RAND_MAX)) < 0.5) {
				return c;
			}
			++c;
		}
	}

	// perform t0 > t1 > t2
	// 3/6 : 2/6 : 1/6 ???
	int select_transition_rule() {
		//return get_random(4);
		double r = rand() / double(RAND_MAX);

        //return 0;

		if(r < 0.5)
			return 0;

		if(r < 0.833)
			return 1;

		if((rand() / double(RAND_MAX)) < 0.5)
			return 2;
		else
			return 3;
	}

	inline unsigned get_random(unsigned i) {
		return rand() % i;
	}

	unsigned get_random_person() {
		return get_random(ped->num_members());
	}

	unsigned get_random_nonleaf() {
		unsigned i = get_random(ped->num_members() - ped->num_leaves());

		while(ped->get_by_index(i)->isleaf())
			++i;

		return i;
	}

	unsigned get_random_nonfounder() {
		return get_random(ped->num_members() - ped->num_founders()) + ped->num_founders();
	}

	unsigned get_random_locus() {
		return get_random(ped->num_markers());
	}

	enum parentage get_random_parent() {
		return static_cast<enum parentage>(get_random(2));
	}

	// select random individual, locus, maternal or paternal, 
	// switch 0 -> 1 or 1 -> 0
	void transition_t0(unsigned locus) {
		transition_t0(get_random_nonfounder(), locus, get_random_parent());
	}

	void transition_t0(unsigned id, unsigned locus, enum parentage p) {
		flip_bit(id, locus, p);
	}

	// select random individual, locus, 
	// for each child, 
	//		if you inherit from maternal change to paternal + vice versa
	void transition_t1(unsigned locus) {		
		transition_t1(ped->get_by_index(get_random_nonleaf()), locus);
	}

	// select random individual, locus, 
	// for each child, 
	//		if you inherit from maternal change to paternal + vice versa		
	void transition_t1(Person* p, unsigned locus) {
		for(unsigned i = 0; i < p->num_children(); ++i) {
			Person* k = p->get_child(i);

			transition_t0(k->get_internalid(), locus, MATERNAL);
			transition_t0(k->get_internalid(), locus, PATERNAL);
		}
	}
	
	// see appendix a1 of simwalk paper for implementing t2a + t2b in terms
	// of t0 and t1
	void transition_t2a(unsigned locus) {
		transition_t2(get_random_nonleaf(), locus, false);
	}

	void transition_t2b(unsigned locus) {
		transition_t2(get_random_nonleaf(), locus, true);
	}

	void transition_t2(unsigned id, unsigned locus, bool same_gender) {
		Person* p = ped->get_by_index(id);
		Person* s = p->num_mates() == 1 ? \
					p->get_mate(0) : \
					p->get_mate(get_random(p->num_mates()));
		
		
		for(unsigned i = 0; i < p->num_children(); ++i) {
			Person* k = p->get_child(i);
			
			// check both parents are correct
			if(s->get_internalid() != (p->ismale() ? k->get_maternalid() : k->get_paternalid()))
				continue;
			
			if(same_gender) {
				if(get(k->get_internalid(), locus, MATERNAL) == \
				   get(k->get_internalid(), locus, PATERNAL)) {
					transition_t0(k->get_internalid(), locus, MATERNAL);
					transition_t0(k->get_internalid(), locus, PATERNAL);
				}
			}
			else {
				if(get(k->get_internalid(), locus, MATERNAL) != \
				   get(k->get_internalid(), locus, PATERNAL)) {
					transition_t0(k->get_internalid(), locus, MATERNAL);
					transition_t0(k->get_internalid(), locus, PATERNAL);
				}
			}
			
			if(not p->isleaf())
			    transition_t1(k, locus);
		}
	}
};

#endif

