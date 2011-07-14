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
        //temperatures.push_back(i / static_cast<double>(num_chains-1));
    }
    
    double tmp = 100.0;
    temperatures.resize(num_chains);
    temperatures[num_chains - 1] = tmp / 100.0;
    for(int i = (num_chains - 2); i > -1; --i) {
        tmp *= 0.95;
        temperatures[i] = tmp / 100.0;
    }
    
    for(unsigned i = 0; i < num_chains; ++i) {
        printf("temperature %d = %f\n", i, temperatures[i]);
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
    /*
    double r = \
        (
            (ls1->likelihood(temp2) + ls2->likelihood(temp1)) - \
            (ls1->likelihood(temp1) + ls2->likelihood(temp2)) \
        );
        
    return log(get_random()) < r;
    */
    
    double prob1 = ls1->likelihood(temp1);
    double prob2 = ls2->likelihood(temp2);
    
    double r = log(get_random());
    double delta_t = (1 / (temp2 * 100)) - (1 / (temp1 * 100));
    
    printf("r = %e, dt = %e, p = %e\n", r, delta_t, prob2-prob1);

    delta_t = log(delta_t);
    
    return r < (prob2 - prob1) * delta_t;
    
}

Peeler* ParallelTempering::run(unsigned iterations) {
    
    unsigned burnin = iterations * 0.1;
    unsigned burst_len = 20;
    
    unsigned total = 0;
    unsigned accepted = 0;
    vector<unsigned> totals(chains.size());
    vector<unsigned> accepteds(chains.size());
    
    for(unsigned i = 0; i < chains.size(); ++i) {
        chains[i]->set_burnin(burnin);
        //chains[i]->anneal(10000);
        
        totals[i] = 0;
        accepteds[i] = 0;
    }
    
    for(unsigned i = 0; i < iterations; i += burst_len) {
    
        for(unsigned j = 0; j < chains.size(); ++j) {
            printf("running chain %d...\n", j);
            
            chains[j]->run(i, burst_len, temperatures[j], peel);
            
            //printf("CHAIN %d %f\n", j, chains[j]->likelihood(0.0));
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
    
    fprintf(stderr, "\n\nfrequency of replica exchanges: %.4f\n\n", accepted / double(total));
    
    for(unsigned i = 1; i < chains.size(); ++i) {
        fprintf(stderr, "chain %d --> %d acceptance freq = %.4f\n", i-1, i, accepteds[i] / double(totals[i]));
    }

    return &peel;
}

