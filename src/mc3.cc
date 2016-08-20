#include <vector>
#include <algorithm>

#include <fstream>
#include <iomanip>

#include "descent_graph.h"
#include "peel_sequence_generator.h"
#include "lod_score.h"
#include "peeler.h"
#include "markov_chain.h"
#include "mc3.h"
#include "progress.h"
#include "sequential_imputation.h"
#include "omp_facade.h"

using namespace std;


void Mc3::_init() {

    lod = new LODscores(map);

    for(int i = 0; i < get_max_threads(); ++i) {
        Peeler* tmp = new Peeler(ped, map, psg, lod, options.sex_linked);
        peelers.push_back(tmp);
    }

    lod->set_trait_prob(peelers[0]->calc_trait_prob());

    for(int i = 0; i < options.mc3_number_of_chains; ++i) {
        double temperature;
        
        if(not options.mc3) {
            temperature = 1.0;
        }
        else if(options.mc3_temperatures.size() > 0) {
            temperature = options.mc3_temperatures[i];
        }
        else {
            temperature = i == 0 ? 1.0 : 1.0 / (1.0 + (0.001 * pow(2.0, i)));
        }

        /*
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
        */

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

LODscores* Mc3::run() {

#define MC3_INFO
#ifdef MC3_INFO
    ofstream f;
    f.open("log");
    f << "iteration chain likelihood coldlikelihood\n";
#endif

    vector<DescentGraph> graphs;
    SequentialImputation si(ped, map, psg, options.sex_linked);

    for(unsigned i = 0; i < chains.size(); ++i) {
        DescentGraph tmp(ped, map, options.sex_linked);

        if(options.si_iterations != 0)
            si.parallel_run(tmp, options.si_iterations);
        else 
            tmp.random_descentgraph();

        graphs.push_back(tmp);
    }


    vector<int> swap_success(chains.size(), 0);
    vector<int> swap_failure(chains.size(), 0);

    // to ensure the progree counter is updated
    if((not options.mc3) or (chains.size() == 1)) {
        options.mc3_exchange_period = 10;
    }

    int spurts = (options.burnin + options.iterations) / options.mc3_exchange_period;

    Progress p((not options.mc3) or (chains.size() == 1) ? "MCMC: " : "MC3: ", spurts);

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

        if(options.mc3) {
            // metropolis step
            int rand = (chains.size() - 1) * get_random();

            double xx = chains[rand]->get_likelihood(graphs[rand]); 
            double yy = chains[rand+1]->get_likelihood(graphs[rand+1]);
            double xy = chains[rand]->get_likelihood(graphs[rand+1]);
            double yx = chains[rand+1]->get_likelihood(graphs[rand]);
            double ratio = (xy + yx) - (xx + yy);
        
            double r = get_random();

            if((r == 0.0) or (log(r) < min(0.0, ratio))) {
                swap(graphs[rand], graphs[rand+1]);
                swap_success[rand] += 1;
            }
            else {
                swap_failure[rand] += 1;
            }
        }
    }

    p.finish();

    // output a file containing the exchange rates
    // between chains
    if(options.mc3) {
        if(options.exchange_filename != "") {
            ofstream ef;
            ef.open(options.exchange_filename.c_str());

            for(int i = 0; i < int(chains.size())-1; ++i) {
                ef << i << " " << swap_success[i] / static_cast<double>(swap_success[i] + swap_failure[i]) << "\n";
            }

            ef.close();
        }

        for(int i = 0; i < int(chains.size())-1; ++i) {
            fprintf(stderr, "%d -- %d : %.3f (%d/%d)\n", \
                i, i+1, swap_success[i] / static_cast<double>(swap_success[i] + swap_failure[i]), \
                swap_success[i], swap_success[i] + swap_failure[i]);
        }
    }

#ifdef MC3_INFO
    f.close();
#endif

    if((int(options.mc3_temperatures.size()) > 0) and (options.mc3_temperatures[0] != 1.0)) {
        exit(EXIT_SUCCESS);
    }

    if((not options.mc3) and (chains.size() != 1)) {
        LODscores* tmp = chains[0]->get_result();
        for(unsigned i = 1; i < chains.size(); ++i) {
            tmp->merge_results(chains[i]->get_result());
        }
        return tmp;
    }

    return chains[0]->get_result();
}

