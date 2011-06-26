#ifndef LKG_SAMPLER_H_
#define LKG_SAMPLER_H_

using namespace std;

#include <cstdlib>


class Pedigree;
class GeneticMap;
class DescentGraph;
class DescentGraphDiff;

class Sampler {

 protected :
    
    Pedigree* ped;
    GeneticMap* map;
    DescentGraph* dg;    
    
    unsigned get_random(unsigned i) {
    	return random() % i;
    }

    unsigned get_random_locus() {
	    return get_random(ped->num_markers());
    }
    
 public :
    Sampler(Pedigree* ped, GeneticMap* map, DescentGraph* dg) : 
        ped(ped), 
        map(map),
        dg(dg) {}
        
    Sampler(const Sampler& rhs) :
        ped(rhs.ped),
        map(rhs.map),
        dg(rhs.dg) {}
    
	virtual ~Sampler() {}
    
    virtual void step(DescentGraphDiff& dgd)=0;
    
    Sampler& operator=(const Sampler& rhs) {
        
        if(&rhs != this) {
            ped = rhs.ped;
            map = rhs.map;
            dg = rhs.dg;
        }
        
        return *this;
    }
};

#endif
