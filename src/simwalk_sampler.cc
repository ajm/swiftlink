using namespace std;

#include <cstdlib>
#include <vector>

#include "genotype.h"
#include "pedigree.h"
#include "person.h"
#include "descent_graph.h"
#include "simwalk_sampler.h"
#include "descent_graph_diff.h"


unsigned SimwalkSampler::_geometric() {
	unsigned c = 1;

	while(1) {
		if((random() / double(RAND_MAX)) < 0.5) {
			return c;
		}
		++c;
	}
}

// go geom(0.5) transitions per step of the chain
// select people + transition rules at random, but not necessarily uniform
// 		perform t0 > t1 > t2
// 		select untyped more than typed
// person + locus need in all rules, add as parameters
	
// t0 - non-founders (founders are a waste of time)
// t1 - non-leaves
// t2 - non-leaves (people with spouses)
void SimwalkSampler::step(DescentGraphDiff& dgd) {
    unsigned steps = _geometric();
	unsigned l = get_random_locus();
	
	dgd.clear();
    
    // corner case : even number of transitions where one transition undoes the other...
    while(dgd.num_transitions() == 0) {
        
        for(unsigned i = 0; i < steps; ++i) {
        
            // could select a neighbouring individual, I don't think sw does this 
            // for the random walk case
            // (a neighbor is a parent, sib, spouse, child or the individual themselves again)
            unsigned p = get_random_person();
            
            //transition_t0(dgd, p,l);
            //return;

            switch(select_transition_rule(p)) {
                case 0 :
                    transition_t0(dgd,p,l);
                    break;

                case 1 :
                    transition_t1(dgd,p,l);
                    break;

                case 2 :
                    transition_t2a(dgd,p,l);
                    break;

                case 3 :
                    transition_t2b(dgd,p,l);
                    break;

                default:
                    break;
            }
        }
    }
}

/*
unsigned SimwalkSampler::select_next_locus(unsigned locus) {
    double r = random() / double(RAND_MAX);
    
    if(r < 0.3333) {
        return locus == 0 ? locus : locus - 1;
    }
    else if(r < 0.6666) {
        return locus;
    }
    else {
        return locus == (ped.num_markers() - 1) ? locus : locus + 1;
    }
}
*/

int SimwalkSampler::select_transition_rule(unsigned person) {
	Person* p = ped.get_by_index(person);
	
	if(p->isleaf()) { // t1 & t2 require spouses
	    type0_count++;
	    return 0;
	}
	
	if(not p->isfounder()) {
        if((random() / double(RAND_MAX)) < 0.75) {
            type0_count++;
            return 0;
        }
    }
    
	if((random() / double(RAND_MAX)) < 0.6) {
	    type1_count++;
	    return 1;
	}	
	
	type2_count++;
	
	if((random() / double(RAND_MAX)) < 0.5) {
		return 2;
	}
	else {
		return 3;
    }
}

unsigned SimwalkSampler::get_random(unsigned i) {
	return random() % i;
}

unsigned SimwalkSampler::get_random_person() {
	return get_random(ped.num_members());
}

unsigned SimwalkSampler::get_random_nonfounder() {
	return get_random(ped.num_members() - ped.num_founders()) + ped.num_founders();
}

unsigned SimwalkSampler::get_random_locus() {
	return get_random(ped.num_markers());
}

enum parentage SimwalkSampler::get_random_parent() {
	return static_cast<enum parentage>(get_random(2));
}

// select random individual, locus, maternal or paternal, 
// switch 0 -> 1 or 1 -> 0
void SimwalkSampler::transition_t0(DescentGraphDiff& dgd, unsigned person, unsigned locus) {
    transition_t0(dgd, person, locus, get_random_parent());
}

void SimwalkSampler::transition_t0(DescentGraphDiff& dgd, unsigned person, unsigned locus, enum parentage parent) {
    dgd.add_transition(person, locus, parent);
}

// select random individual, locus, 
// for each child, 
//		if you inherit from maternal change to paternal + vice versa		
void SimwalkSampler::transition_t1(DescentGraphDiff& dgd, unsigned person, unsigned locus) {
    Person* p = ped.get_by_index(person);

	for(unsigned i = 0; i < p->num_children(); ++i) {
		Person* k = p->get_child(i);
		
		transition_t0(dgd, k->get_internalid(), locus, MATERNAL);
		transition_t0(dgd, k->get_internalid(), locus, PATERNAL);
	}
}
	
// see appendix a1 of simwalk paper for implementing t2a + t2b in terms
// of t0 and t1
void SimwalkSampler::transition_t2a(DescentGraphDiff& dgd, unsigned person, unsigned locus) {
	transition_t2(dgd, person, locus, false);
}

void SimwalkSampler::transition_t2b(DescentGraphDiff& dgd, unsigned person, unsigned locus) {
	transition_t2(dgd, person, locus, true);
}

void SimwalkSampler::transition_t2(DescentGraphDiff& dgd, unsigned id, unsigned locus, bool same_gender) {
	Person* p = ped.get_by_index(id);
	Person* s = p->num_mates() == 1 ? \
				p->get_mate(0) : \
				p->get_mate(get_random(p->num_mates()));
	
	
	for(unsigned i = 0; i < p->num_children(); ++i) {
		Person* k = p->get_child(i);
		
		// ignore step-children
		if(s->get_internalid() != (p->ismale() ? k->get_maternalid() : k->get_paternalid()))
			continue;
		
		if(same_gender) {
			if(dg.get(k->get_internalid(), locus, MATERNAL) == \
			   dg.get(k->get_internalid(), locus, PATERNAL)) {
			   
				transition_t0(dgd, k->get_internalid(), locus, MATERNAL);
				transition_t0(dgd, k->get_internalid(), locus, PATERNAL);
			}
		}
		else {
			if(dg.get(k->get_internalid(), locus, MATERNAL) != \
			   dg.get(k->get_internalid(), locus, PATERNAL)) {
			   
				transition_t0(dgd, k->get_internalid(), locus, MATERNAL);
				transition_t0(dgd, k->get_internalid(), locus, PATERNAL);
			}
		}
		
		if(not p->isleaf())
		    transition_t1(dgd, k->get_internalid(), locus);
	}
}

