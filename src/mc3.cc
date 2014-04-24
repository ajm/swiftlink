using namespace std;

#include <vector>
#include <algorithm>

#include <fstream>
#include <iomanip>

#include "omp.h"

#include "descent_graph.h"
#include "peel_sequence_generator.h"
#include "lod_score.h"
#include "peeler.h"
#include "markov_chain.h"
#include "mc3.h"
#include "progress.h"
#include "sequential_imputation.h"

void Mc3::_init() {

    lod = new LODscores(map);

    for(int i = 0; i < omp_get_max_threads(); ++i) {
        Peeler* tmp = new Peeler(ped, map, psg, lod);
        peelers.push_back(tmp);
    }

    lod->set_trait_prob(peelers[0]->calc_trait_prob());

    for(int i = 0; i < options.mc3_number_of_chains; ++i) {
        //double temperature = 1.0 / (1.0 + (options.mc3_temperature * i));
        double temperature = 1.0;
        //if(i > 0)
        //    temperature = 1.0 / (1.0 + (options.mc3_temperature * pow(2.0, i)));

        switch(i) {
            case 0 : temperature = 1.000; break;
            case 1 : temperature = 0.998; break;
            case 2 : temperature = 0.996; break;
            case 3 : temperature = 0.992; break;
            case 4 : temperature = 0.984; break;
            case 5 : temperature = 0.969; break;
            case 6 : temperature = 0.940; break;
            case 7 : temperature = 0.887; break;
            case 8 : temperature = 0.820; break;
            case 9 : temperature = 0.760; break;
            case 10: temperature = 0.700; break;
            case 11: temperature = 0.640; break;
            case 12: temperature = 0.580; break;
            case 13: temperature = 0.520; break;
            case 14: temperature = 0.460; break;
            case 15: temperature = 0.400; break;
            default:
                abort();
        }


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

#define MC3_INFO
#ifdef MC3_INFO
    ofstream f;
    f.open("log");
    f << "iteration chain likelihood coldlikelihood\n";
#endif

    if(chains.size() == 1) {
        Progress p("MCMC: ", (options.burnin + options.iterations) / 10);

        for(int i = 0; i < options.burnin + options.iterations; i += 10) {
            chains[0]->step(dg, i, 10);

            p.increment();

#ifdef MC3_INFO
            double lik = chains[0]->get_likelihood(dg);
            f << i << " 0 " << lik << " " << lik << "\n";
#endif
        }

        p.finish();

#ifdef MC3_INFO
        f.close();
#endif

        return chains[0]->get_result();
    }


    vector<int> swap_success(chains.size(), 0);
    vector<int> swap_failure(chains.size(), 0);
/*
    if(dg.get_likelihood() == LOG_ILLEGAL) {
        fprintf(stderr, "error: descent graph illegal pre-markov chain...\n");
        abort();
    }
*/

    vector<DescentGraph> graphs;
    SequentialImputation si(ped, map, psg);
    
    for(unsigned i = 0; i < chains.size(); ++i) {
        DescentGraph tmp(ped, map);
        
        si.parallel_run(tmp, 1000);

        graphs.push_back(tmp);
    }

    int spurts = (options.burnin + options.iterations) / options.mc3_exchange_period;

    Progress p("MC3: ", spurts);

    for(int i = 0; i < spurts; ++i) {
        // advance all chains
        for(unsigned j = 0; j < chains.size(); ++j) {
            chains[j]->step(graphs[j], i * spurts, options.mc3_exchange_period);

#ifdef MC3_INFO
            f << (i * options.mc3_exchange_period) \
              << " " \
              << j \
              << " " \
              << chains[j]->get_likelihood(graphs[j]) \
              << " " \
              << chains[0]->get_likelihood(graphs[j]) \
              << "\n";
#endif
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
        fprintf(stderr, "%d -- %d : %.3f (%d/%d)\n", \
                i, i+1, swap_success[i] / static_cast<double>(swap_success[i] + swap_failure[i]), \
                swap_success[i], swap_success[i] + swap_failure[i]);
    }

#ifdef MC3_INFO
    f.close();
#endif

    return chains[0]->get_result();
}

