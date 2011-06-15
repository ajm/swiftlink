#ifndef LKG_DESCENTGRAPHDIFF_H_
#define LKG_DESCENTGRAPHDIFF_H_

#include <vector>
#include <algorithm>
#include <sstream>

#include "descent_graph_types.h"


class Transition {
    
    unsigned person;
    unsigned locus;
    enum parentage parent;
        
 public:
    Transition(unsigned person, unsigned locus, enum parentage parent) 
        : person(person), locus(locus), parent(parent) {}
    
    unsigned get_person() { return person; }
    unsigned get_locus() { return locus; }
    enum parentage get_parent() { return parent; }
    
    bool operator==(const Transition& t) const {
		return (person == t.person) and \
		       (locus == t.locus) and \
		       (parent == t.parent);
	}
	
	string debug_string() {
	    stringstream ss;
	    
	    ss << "Person: " << person << ", ";
	    ss << "Locus: " << locus << ", ";
	    ss << "Allele: " << parent;
	    
	    return ss.str();
	}
};

// XXX to start off I am going to assume that all 
// diffs are 1 step t0s and can therefore be described
// as a person,locus,parent tuple
class DescentGraphDiff {
    
    vector<Transition> transitions;
    double sum_prior_prob;
    double prob;
    int recombinations;

 public :        
    DescentGraphDiff() :
        transitions(), sum_prior_prob(0.0), prob(0.0), recombinations(-1) {}        
    ~DescentGraphDiff() {}
    
    void clear() {
        transitions.clear();
    }
    
    void add_transition(unsigned person, unsigned locus, enum parentage parent) {
        Transition t(person, locus, parent);
        vector<Transition>::iterator it;
        
        // if transition already exists, then transitioning again will undo it
        it = find(transitions.begin(), transitions.end(), t);
        
        if(it == transitions.end()) {
            transitions.push_back(t);
        }
        else {
            transitions.erase(it);
        }
    }
    
    unsigned num_transitions() { return transitions.size(); }
    Transition& get_transition(unsigned i) { return transitions[i]; }
    
    double get_sumprior() { return sum_prior_prob; }
    void set_sumprior(double p) { sum_prior_prob = p; }
    
    double get_prob() { return prob; }
    void set_prob(double p) { prob = p; }
    
    int get_recombinations() { return recombinations; }
    void set_recombinations(int r) { recombinations = r; }
    
    void print() {
        for(unsigned i = 0; i < transitions.size(); ++i) {
            printf("%s\n", transitions[i].debug_string().c_str());
        }
    }
};

#endif

