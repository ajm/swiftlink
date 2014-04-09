using namespace std;

#include <vector>
#include <algorithm>

#include "omp.h"

#include "descent_graph.h"
#include "peel_sequence_generator.h"
#include "lod_score.h"
#include "peeler.h"
#include "markov_chain.h"
#include "mc3.h"
#include "progress.h"

void Mc3::_init() {

    lod = new LODscores(map);

    for(int i = 0; i < omp_get_max_threads(); ++i) {
        Peeler* tmp = new Peeler(ped, map, psg, lod);
        peelers.push_back(tmp);
    }

    lod->set_trait_prob(peelers[0]->calc_trait_prob());

    for(int i = 0; i < options.mc3_number_of_chains; ++i) {
        double temperature = 1.0 / (1.0 + (options.mc3_temperature * i));
        
        fprintf(stderr, "Creating Markov chain %d, temperature = %.3f\n", i, temperature);

        MarkovChain* tmp = new MarkovChain(ped, map, psg, options, temperature);
        chains.push_back(tmp);
    }
}

void Mc3::_kill() {
    for(int i = 0; i < int(peelers.size()); ++i) {
        delete peelers[i];
    }
    for(int i = 0; i < int(chains.size()); ++i) {
        delete chains[i];
    }
}

LODscores* Mc3::run(DescentGraph& dg) {

    vector<int> swap_success(chains.size(), 0);
    vector<int> swap_failure(chains.size(), 0);

    if(dg.get_likelihood() == LOG_ILLEGAL) {
        fprintf(stderr, "error: descent graph illegal pre-markov chain...\n");
        abort();
    }

    // make a copy of the starting state for each chain
    vector<DescentGraph> graphs;
    for(unsigned i = 0; i < chains.size(); ++i) {
        graphs.push_back(dg);
    }

    int spurts = (options.burnin + options.iterations) / options.mc3_exchange_period;

    Progress p("MC3: ", spurts);

    for(int i = 0; i < spurts; ++i) {
        // advance all chains
        for(unsigned j = 0; j < chains.size(); ++j) {
            chains[j]->step(graphs[j], i * spurts, options.mc3_exchange_period);
        }

        p.increment();

        // metropolis step
        int rand = (chains.size() - 1) * get_random();

        double xx = chains[rand]->get_likelihood(graphs[rand]); 
        double yy = chains[rand+1]->get_likelihood(graphs[rand+1]);
        double xy = chains[rand]->get_likelihood(graphs[rand+1]);
        double yx = chains[rand+1]->get_likelihood(graphs[rand]);

        double r = get_random();

        if((r == 0.0) or (log(r) < min(0.0, (xy + yx) - (xx + yy)))) {
            swap(graphs[rand], graphs[rand+1]);
            swap_success[rand] += 1;
        }
        else {
            swap_failure[rand] += 1;
        }
    }

    p.finish();

    for(int i = 0; i < int(chains.size())-1; ++i) {
        fprintf(stderr, "%d -- %d : %.3f\n", i, i+1, swap_success[i] / static_cast<double>(swap_success[i] + swap_failure[i]));
    }

    return chains[0]->get_result();
}

