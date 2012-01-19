#ifndef LKG_LOCUSSAMPLER_H_
#define LKG_LOCUSSAMPLER_H_

using namespace std;

#include <vector>

#include "descent_graph.h"
#include "trait.h"
#include "sampler.h"
#include "sampler_rfunction.h"


class Pedigree;
class GeneticMap;
class PeelSequenceGenerator;


class LocusSampler : Sampler {
    
    vector<SamplerRfunction> rfunctions;
    unsigned int locus;
    
    void init_rfunctions(PeelSequenceGenerator& psg);    
    unsigned sample_mi(DescentGraph& dg, enum trait allele, enum phased_trait trait, unsigned personid, enum parentage parent);
    unsigned sample_homo_mi(DescentGraph& dg, unsigned personid, enum parentage parent);
    unsigned sample_hetero_mi(enum trait allele, enum phased_trait trait);
    
    void sample_meiosis_indicators(vector<int>& pmk, DescentGraph& dg);
    
    void set_all(bool left, bool right);

    
 public :
    LocusSampler(Pedigree* ped, GeneticMap* map, PeelSequenceGenerator& psg, unsigned int locus) :
        Sampler(ped, map), 
        rfunctions(),
        locus(locus) {
        
        init_rfunctions(psg);
    }
    
    ~LocusSampler() {}
    
    LocusSampler(const LocusSampler& rhs) : 
        Sampler(rhs.ped, rhs.map),
        rfunctions(rhs.rfunctions),
        locus(rhs.locus) {}
    
    LocusSampler& operator=(const LocusSampler& rhs) {
        
        if(this != &rhs) {
            Sampler::operator=(rhs);
            rfunctions = rhs.rfunctions;
            locus = rhs.locus;
        }
        
        return *this;        
    }
    
    virtual void step(DescentGraph& dg, unsigned parameter);
    
    void sequential_imputation(DescentGraph& dg);
    void reset();
    
    void set_locus(unsigned int locus) {
    
        this->locus = locus;
    
        for(unsigned int i = 0; i < rfunctions.size(); ++i) {
            rfunctions[i].set_locus(locus);
        }
    }
};

#endif

