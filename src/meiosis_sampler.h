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
    
    FounderAlleleGraph  f1;
    FounderAlleleGraph2 f2;
    FounderAlleleGraph3 f3;
    FounderAlleleGraph4 f4;
    MeiosisMatrix* matrix;
    vector<int> seq;
    
    void init_matrices();
    void copy_matrices(const MeiosisSampler& rhs);
    void kill_matrices();
    double graph_likelihood(DescentGraph& dg, unsigned person_id, unsigned locus, enum parentage parent, unsigned value);
    double initial_likelihood(DescentGraph& dg, unsigned locus);
    void incremental_likelihood(DescentGraph& dg, unsigned person_id, unsigned locus, enum parentage parent, double* meiosis0, double* meiosis1);
    void find_founderallelegraph_ordering();
    
 public :
    MeiosisSampler(Pedigree* ped, GeneticMap* map) :
        Sampler(ped, map), 
        f1(ped, map),
        f2(ped, map, 0),
        f3(ped, map),
        f4(ped, map, 0),
        matrix(NULL),
        seq() {
    
        init_matrices();
        
        find_founderallelegraph_ordering();
        f3.set_sequence(&seq);
        f4.set_sequence(&seq);
    }
    
    MeiosisSampler(const MeiosisSampler& rhs) :
        Sampler(rhs.ped, rhs.map),
        f1(rhs.f1),
        f2(rhs.f2),
        f3(rhs.f3),
        f4(rhs.f4),
        matrix(NULL),
        seq(rhs.seq) {
        
        init_matrices();
        copy_matrices(rhs);
    }
    
    virtual ~MeiosisSampler() {
        kill_matrices();
    }
    
    MeiosisSampler& operator=(const MeiosisSampler& rhs) {
        
        if(this != &rhs) {
            Sampler::operator=(rhs);
            copy_matrices(rhs);
            seq = rhs.seq;
        }
        
        return *this;
    }
    
    virtual void step(DescentGraph& dg, unsigned parameter);
};

#endif

