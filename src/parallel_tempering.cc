using namespace std;

#include <vector>
#include <cstdlib>
#include <cmath>

#include "parallel_tempering.h"
#include "locus_sampler.h"


void ParallelTempering::init_chains(unsigned num_chains) {
    //chains.resize(num_chains);
    //temperatures.resize(num_chains);
    
    for(unsigned i = 0; i < num_chains; ++i) {
        LocusSampler* s = new LocusSampler(ped, map);
        chains.push_back(s);
        temperatures.push_back(i / static_cast<double>(num_chains));
    }
}

void ParallelTempering::copy_chains(const ParallelTempering& rhs) {
    chains.clear();
    for(unsigned i = 0; i < rhs.chains.size(); ++i) {
        LocusSampler* s = new LocusSampler(*(rhs.chains[i]));
        chains.push_back(s);
    }
}

void ParallelTempering::kill_chains() {
    for(unsigned i = 0; i < chains.size(); ++i) {
        delete chains[i];
    }
}

double ParallelTempering::get_random() {
    return random() / static_cast<double>(RAND_MAX);
}

bool ParallelTempering::exchange_replicas(LocusSampler* ls1, LocusSampler* ls2, double temp1, double temp2) {
    double r = \
        (
            (ls1->likelihood(temp2) + ls2->likelihood(temp1)) - \
            (ls1->likelihood(temp1) + ls2->likelihood(temp2)) \
        );
/*
    fprintf(stderr, "(%.4f)(%.4f)/(%.4f)(%.4f) = %f  log(%f)\n", 
        ls1->likelihood(temp2), 
        ls2->likelihood(temp1), 
        ls1->likelihood(temp1), 
        ls2->likelihood(temp2), 
        r, exp(r));
*/
    return log(get_random()) < r;
}

Peeler* ParallelTempering::run(unsigned iterations) {
    
    unsigned burnin = iterations * 0.1;
    unsigned burst_len = 5;
    
    unsigned total = 0;
    unsigned accepted = 0;
    vector<unsigned> totals(chains.size());
    vector<unsigned> accepteds(chains.size());
    
    for(unsigned i = 0; i < chains.size(); ++i) {
        chains[i]->set_burnin(burnin);
        chains[i]->anneal(10000);
        
        totals[i] = 0;
        accepteds[i] = 0;
    }
    
    for(unsigned i = 0; i < iterations; i += burst_len) {
    
        for(unsigned j = 0; j < chains.size(); ++j) {
            printf("running chain %d...\n", j);
        
            chains[j]->run(i, burst_len, temperatures[j], peel);
        }
        
        //for(unsigned j = 1; j < chains.size(); ++j) {
        for(unsigned j = chains.size() - 1; j > 0; --j) {
            total++;
            totals[j] += 1;
            
            if(exchange_replicas(chains[j-1], chains[j], temperatures[j-1], temperatures[j])) {
                printf("\texchanged %d (%.4f) and %d (%.4f)\n", j-1, temperatures[j-1], j, temperatures[j]);
                accepted++;
                accepteds[j] += 1;
            
                LocusSampler* tmp = chains[j-1];
                chains[j-1] = chains[j];
                chains[j] = tmp;
            }
        }
    }
    
    printf("\n\nfrequency of replica exchanges: %.4f\n\n", accepted / double(total));
    
    for(unsigned i = 1; i < chains.size(); ++i) {
        printf("chain %d --> %d acceptance freq = %.4f\n", i-1, i, accepteds[i] / double(totals[i]));
    }
    
    return &peel;
}

