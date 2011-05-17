#ifndef LKG_DESCENTGRAPHDIFF_H_
#define LKG_DESCENTGRAPHDIFF_H_

#include <vector>

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
    DescentGraphDiff() {}        
    ~DescentGraphDiff() {}
    
    void add_transition(unsigned person, unsigned locus, enum parentage parent) {
        Transition t(person, locus, parent);
        transitions.push_back(t);
    }
    
    unsigned num_transitions() { return transitions.size(); }
    Transition& get_transition(unsigned i) { return transitions[i]; }
    
    double get_sumprior() { return sum_prior_prob; }
    void set_sumprior(double p) { sum_prior_prob = p; }
    
    double get_prob() { return prob; }
    void set_prob(double p) { prob = p; }
    
    int get_recombinations() { return recombinations; }
    void set_recombinations(int r) { recombinations = r; }
};

#endif

