#ifndef LKG_LOCUSSAMPLER_H_
#define LKG_LOCUSSAMPLER_H_

#include <vector>

#include "descent_graph.h"
#include "trait.h"
#include "sampler.h"
#include "sampler_rfunction.h"

using namespace std;


class Pedigree;
class GeneticMap;
class PeelSequenceGenerator;


class LocusSampler : Sampler {
    
    vector<SamplerRfunction> rfunctions;
    unsigned int locus;
    bool ignore_left;
    bool ignore_right;
    bool sex_linked;
    
    void init_rfunctions(PeelSequenceGenerator* psg);    
    unsigned sample_mi(DescentGraph& dg, enum trait allele, enum phased_trait trait, unsigned personid, enum parentage parent);
    unsigned sample_homo_mi(DescentGraph& dg, unsigned personid, enum parentage parent);
    unsigned sample_hetero_mi(enum trait allele, enum phased_trait trait);
    
    void sample_meiosis_indicators(vector<int>& pmk, DescentGraph& dg);

    
 public :
    LocusSampler(Pedigree* ped, GeneticMap* map, PeelSequenceGenerator* psg, unsigned int locus, bool sex_linked) :
        Sampler(ped, map), 
        rfunctions(),
        locus(locus),
        ignore_left(false),
        ignore_right(false),
        sex_linked(sex_linked) {
        
        init_rfunctions(psg);
    }
    
    ~LocusSampler() {}
    
    LocusSampler(const LocusSampler& rhs) : 
        Sampler(rhs.ped, rhs.map),
        rfunctions(rhs.rfunctions),
        locus(rhs.locus),
        ignore_left(rhs.ignore_left),
        ignore_right(rhs.ignore_right) {}
    
    LocusSampler& operator=(const LocusSampler& rhs) {
        
        if(this != &rhs) {
            Sampler::operator=(rhs);
            rfunctions = rhs.rfunctions;
            locus = rhs.locus;
            ignore_left = rhs.ignore_left;
            ignore_right = rhs.ignore_right;
        }
        
        return *this;        
    }
    
    virtual void step(DescentGraph& dg, unsigned parameter);
    
    double locus_by_locus(DescentGraph& dg);
    double sequential_imputation(DescentGraph& dg);
    double start_from(DescentGraph& dg, unsigned int starting_locus);
    void reset();
    
    void set_locus(unsigned int locus, bool ignore_left, bool ignore_right);
    void set_locus_minimal(unsigned int locus);
};

#endif

