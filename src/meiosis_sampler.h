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
#include "founder_allele_graph2.h"
#include "founder_allele_graph.h"

class Pedigree;
class GeneticMap;


// we are using log likelihoods here
// unlike in the L-sampler
class MeiosisMatrix {
    double prob[2];
    
 public :
    double& operator[](unsigned int i) {
        return prob[i];
    }
    
    const double& operator[](unsigned int i) const {
        return prob[i];
    }
    
    void print() {
        printf("0 %e\n1 %e\n\n", prob[0], prob[1]);
    }
    
    void normalise() {
    
        if((prob[0] == LOG_ZERO) and (prob[1] == LOG_ZERO)) {
            fprintf(stderr, \
                "error: both meiosis likelihoods were zero (%s:%d)\n", \
                __FILE__, __LINE__);
            abort();
        }
    
        if(prob[0] == LOG_ZERO) {
            prob[1] = 0.0;
            return;
        }
        
        if(prob[1] == LOG_ZERO) {
            prob[0] = 0.0;
            return;
        }
        
        double tot = log_sum(prob[0], prob[1]);
        
        prob[0] -= tot;
        prob[1] -= tot;
    }
    
    unsigned sample() {
        return (log(random() / static_cast<double>(RAND_MAX)) < prob[0]) ? 0 : 1;
    }
};

class MeiosisSampler : Sampler {
    
    FounderAlleleGraph3 f3;
    //FounderAlleleGraph4 f4;
    vector<FounderAlleleGraph4> f4;
    vector<double> matrix; // have this double the length
    vector<int> seq;
    
    /*
    void init_matrices();
    void copy_matrices(const MeiosisSampler& rhs);
    void kill_matrices();
    */
    
    double graph_likelihood(DescentGraph& dg, unsigned person_id, unsigned locus, enum parentage parent, unsigned value);
    double initial_likelihood(DescentGraph& dg, unsigned locus);
    void incremental_likelihood(DescentGraph& dg, unsigned person_id, unsigned locus, enum parentage parent, double* meiosis0, double* meiosis1);
    void find_founderallelegraph_ordering();
    
    void normalise(int locus);
    unsigned sample(int locus);
    
 public :
    MeiosisSampler(Pedigree* ped, GeneticMap* map) :
        Sampler(ped, map),
        f3(ped, map),
        f4(map->num_markers(), FounderAlleleGraph4(ped, map)),
        matrix(map->num_markers() * 2),
        seq() {
    
        find_founderallelegraph_ordering();
        
        f3.set_sequence(&seq);
        
        for(unsigned int i = 0; i < map->num_markers(); ++i) {
            f4[i].set_sequence(&seq);
            f4[i].set_locus(i);
        }
    }
    
    MeiosisSampler(const MeiosisSampler& rhs) :
        Sampler(rhs.ped, rhs.map),
        f3(rhs.f3),
        f4(rhs.f4),
        matrix(rhs.matrix),
        seq(rhs.seq) {}
    
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

