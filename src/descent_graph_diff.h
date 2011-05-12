#ifndef LKG_DESCENTGRAPHDIFF_H_
#define LKG_DESCENTGRAPHDIFF_H_

#include "descent_graph.h"

// XXX to start off I am going to assume that all 
// diffs are 1 step t0s and can therefore be described
// as a person,locus,parent tuple
class DescentGraphDiff {

    unsigned person;
    unsigned locus;
    enum parentage parent;

  public :        
    DescentGraphDiff(unsigned person, unsigned locus, enum parentage parent)
        : person(person), locus(locus), parent(parent) {}
        
    ~DescentGraphDiff() {}
    
    unsigned get_person() { return person; }
    unsigned get_locus() { return locus; }
    enum parentage get_parent() { return parent; }
};

#endif

