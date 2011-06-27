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
    
    unsigned sample_mi(unsigned allele, enum phased_trait trait, 
                       unsigned personid, unsigned locus, enum parentage parent);
    unsigned sample_homo_mi(unsigned personid, unsigned locus, enum parentage parent);
    unsigned sample_hetero_mi(unsigned allele, enum phased_trait trait);
    
 public :
    LocusSampler(Pedigree* ped, GeneticMap* map, DescentGraph* dg);        
    
    ~LocusSampler();
    
    LocusSampler(const LocusSampler& rhs) :
        Sampler(rhs), 
        rfunctions(rhs.rfunctions) {}
        
    LocusSampler& operator=(const LocusSampler& rhs) {
        if(this != &rhs) {
			Sampler::operator=(rhs);
			//rfunctions = rhs.rfunctions;
		}        
		return *this;
    }
    
    void step(DescentGraphDiff& dgd);
};

#endif

