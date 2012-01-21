#ifndef LKG_SAMPLER_H_
#define LKG_SAMPLER_H_

using namespace std;

#include <cstdlib>

#include "types.h"
#include "descent_graph.h"
#include "pedigree.h"
#include "genetic_map.h"
#include "random.h"


class Sampler {

 protected :
    Pedigree* ped;
    GeneticMap* map;
    /*
    double get_random() {
        return random() / static_cast<double>(RAND_MAX);
    }
    
    unsigned get_random(unsigned i) {
    	return random() % i;
    }
    */
    unsigned get_random_locus() {
	    return get_random(map->num_markers());
    }
    
    unsigned get_random_nonfounder() {
        return get_random(ped->num_members() - ped->num_founders()) + ped->num_founders();
    }
    
    enum parentage get_random_meiosis() {
        return static_cast<enum parentage>(get_random(2));
    }

 public :
    Sampler(Pedigree* ped, GeneticMap* map) : 
        ped(ped), 
        map(map) {}
        
    Sampler(const Sampler& rhs) :
        ped(rhs.ped),
        map(rhs.map) {}
    
	virtual ~Sampler() {}
    
    Sampler& operator=(const Sampler& rhs) {
        
        if(&rhs != this) {
            ped = rhs.ped;
            map = rhs.map;
        }
        
        return *this;
    }
    
    virtual void step(DescentGraph& dg, unsigned parameter)=0;
};

#endif

