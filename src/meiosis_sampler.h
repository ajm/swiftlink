#ifndef LKG_MEIOSISSAMPLER_H_
#define LKG_MEIOSISSAMPLER_H_

using namespace std;

#include <limits>
#include <vector>

#include "types.h"
#include "sampler.h"
#include "logarithms.h"
#include "founder_allele_graph4.h"

#include "founder_allele_graph3.h"


class Pedigree;
class GeneticMap;


class MeiosisSampler : Sampler {
    
    vector<FounderAlleleGraph4> f4;
    vector<double> matrix;
    vector<int> seq;
    
    FounderAlleleGraph3 f3;
    
    
    double graph_likelihood(DescentGraph& dg, unsigned person_id, unsigned locus, enum parentage parent, unsigned value);
    double initial_likelihood(DescentGraph& dg, unsigned locus);
    void incremental_likelihood(DescentGraph& dg, unsigned person_id, unsigned locus, enum parentage parent, double* meiosis0, double* meiosis1);
    void find_founderallelegraph_ordering();
    
    unsigned sample(int locus);
    
 public :
    MeiosisSampler(Pedigree* ped, GeneticMap* map) :
        Sampler(ped, map),
        f4(map->num_markers(), FounderAlleleGraph4(ped, map)),
        matrix(map->num_markers() * 2),
        seq(),
        f3(ped, map) {
        
        find_founderallelegraph_ordering();
        
        for(unsigned int i = 0; i < map->num_markers(); ++i) {
            f4[i].set_sequence(&seq);
            f4[i].set_locus(i);
        }
        
        f3.set_sequence(&seq);
    }
    
    MeiosisSampler(const MeiosisSampler& rhs) :
        Sampler(rhs.ped, rhs.map),
        f4(rhs.f4),
        matrix(rhs.matrix),
        seq(rhs.seq),
        f3(rhs.f3) {}
    
    virtual ~MeiosisSampler() {}
    
    MeiosisSampler& operator=(const MeiosisSampler& rhs) {
        
        if(this != &rhs) {
            Sampler::operator=(rhs);
            f4 = rhs.f4;
            matrix = rhs.matrix;
            seq = rhs.seq;
        }
        
        return *this;
    }
    
    void reset(DescentGraph& dg);
    
    virtual void step(DescentGraph& dg, unsigned parameter);
};

#endif

