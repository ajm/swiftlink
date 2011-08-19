#ifndef LKG_MEIOSISSAMPLER_H_
#define LKG_MEIOSISSAMPLER_H_

using namespace std;

#include <limits>

#include "sampler.h"
#include "logarithms.h"
#include "descent_graph_types.h"
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
    
    FounderAlleleGraph2 f;
    FounderAlleleGraph fag;
    MeiosisMatrix* matrix;
    
    void init_matrices();
    void copy_matrices(const MeiosisSampler& rhs);
    void kill_matrices();
    double graph_likelihood(DescentGraph& dg, unsigned person_id, unsigned locus, enum parentage parent, unsigned value);
    double initial_likelihood(DescentGraph& dg, unsigned locus);
    void incremental_likelihood(DescentGraph& dg, unsigned person_id, unsigned locus, enum parentage parent, double* meiosis0, double* meiosis1);
    
 public :
    MeiosisSampler(Pedigree* ped, GeneticMap* map) :
        Sampler(ped, map), 
        f(ped, map, 0),
        fag(map, ped),
        matrix(NULL) {
    
        init_matrices();
    }
    
    MeiosisSampler(const MeiosisSampler& rhs) :
        Sampler(rhs.ped, rhs.map),
        f(ped, map, 0),
        fag(map, ped),
        matrix(NULL) {
        
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
        }
        
        return *this;
    }
    
    virtual void step(DescentGraph& dg, unsigned parameter);
};

#endif
