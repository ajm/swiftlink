#ifndef LKG_LOCUSSAMPLER_H_
#define LKG_LOCUSSAMPLER_H_

using namespace std;

#include <vector>

#include "sampler.h"
#include "descent_graph_diff.h"

class Pedigree;
class DescentGraph;
class Rfunction;

class LocusSampler : public Sampler {

    vector<Rfunction*> rfunctions;
    
    LocusSampler(Pedigree* ped, DescentGraph* dg) :
        Sampler(ped, dg) {}
        
    ~LocusSampler() {}
    
    LocusSampler(const LocusSampler& rhs) :
        Sampler(rhs) {}
        
    LocusSampler& operator=(const LocusSampler& rhs) {
        if(this != &rhs) {
			Sampler::operator=(rhs);
		}        
		return *this;
    }
    
    void step(DescentGraphDiff& dgd);
};

#endif

