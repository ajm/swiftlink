#ifndef LKG_LOCUSSAMPLER_H_
#define LKG_LOCUSSAMPLER_H_

using namespace std;

#include <vector>

#include "descent_graph.h"
#include "trait.h"
#include "sampler.h"


class Pedigree;
class GeneticMap;
class SamplerRfunction;
class PeelSequenceGenerator;


class LocusSampler : Sampler {
    
    vector<SamplerRfunction*> rfunctions;
    
    void init_rfunctions(PeelSequenceGenerator& psg);
    void copy_rfunctions(const LocusSampler& rhs);
    void kill_rfunctions();
    
    unsigned sample_mi(DescentGraph& dg, enum trait allele, enum phased_trait trait, unsigned personid, unsigned locus, enum parentage parent);
    unsigned sample_homo_mi(DescentGraph& dg, unsigned personid, unsigned locus, enum parentage parent);
    unsigned sample_hetero_mi(enum trait allele, enum phased_trait trait);
    
    void sample_meiosis_indicators(vector<int>& pmk, DescentGraph& dg, unsigned locus);
    
    void set_all(bool left, bool right);

    
 public :
    LocusSampler(Pedigree* ped, GeneticMap* map, PeelSequenceGenerator& psg) :
        Sampler(ped, map), 
        rfunctions() {
        
        init_rfunctions(psg);
    }
    
    ~LocusSampler() {
        kill_rfunctions();
    }
    
    LocusSampler(const LocusSampler& rhs) : 
        Sampler(rhs.ped, rhs.map),
        rfunctions() {
        
        copy_rfunctions(rhs);
    }
    
    LocusSampler& operator=(const LocusSampler& rhs) {
        
        if(this != &rhs) {
            Sampler::operator=(rhs);
            
            kill_rfunctions();
            copy_rfunctions(rhs);
        }
        
        return *this;        
    }
    
    virtual void step(DescentGraph& dg, unsigned parameter);
    void sequential_imputation(DescentGraph& dg);
    void reset();
};

#endif

