#ifndef LKG_SAMPLER_H_
#define LKG_SAMPLER_H_

using namespace std;

#include <cstdlib>

class Pedigree;
class DescentGraph;
class DescentGraphDiff;

class Sampler {
    
    unsigned get_random(unsigned i) {
    	return random() % i;
    }

    unsigned get_random_locus() {
	    return get_random(ped->num_markers());
    }
    
 protected :
    Pedigree* ped;
    DescentGraph* dg;

 public :
    Sampler(Pedigree* ped, DescentGraph* dg) : 
        ped(ped), 
        dg(dg) {}
        
    Sampler(const Sampler& rhs) :
        ped(rhs.ped),
        dg(rhs.dg) {}
    
	virtual ~Sampler() {}
    virtual void step(DescentGraphDiff& dgd)=0;
    
    Sampler& operator=(const Sampler& rhs) {
        
        if(&rhs != this) {
            ped = rhs.ped;
            dg = rhs.dg;
        }
        
        return *this;
    }
};

#endif

