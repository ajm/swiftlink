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

    vector<SamplerRfunction*> rfunctions;
    
 public :
    LocusSampler(Pedigree* ped, GeneticMap* map, DescentGraph* dg);        
    
    ~LocusSampler() {}
    
    LocusSampler(const LocusSampler& rhs) :
        Sampler(rhs), 
        rfunctions(rhs.rfunctions) {}
        
    LocusSampler& operator=(const LocusSampler& rhs) {
        if(this != &rhs) {
			Sampler::operator=(rhs);
		}        
		return *this;
    }
    
    void step(DescentGraphDiff& dgd);
};

#endif

