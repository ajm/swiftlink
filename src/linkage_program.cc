using namespace std;

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "peel_sequence_generator.h"
#include "simwalk_descent_graph.h"
#include "simulated_annealing.h"
#include "linkage_program.h"
#include "disease_model.h"
#include "markov_chain.h"
#include "genetic_map.h"
#include "peel_matrix.h"
#include "rfunction.h"
#include "pedigree.h"
#include "progress.h"
#include "peeling.h"
#include "peeler.h"


bool LinkageProgram::run() {
    bool ret = true;
    
    fprintf(stderr, "%s\n", dm.debug_string().c_str());
    fprintf(stderr, "%s\n", map.debug_string().c_str());
    
    

    // XXX need to know how to do this properly, 
    // look up better random numbers for simulations etc
    srandom(time(NULL));

    for(unsigned int i = 0; i < pedigrees.size(); ++i) {
        pedigrees[i].print();
        ret &= run_pedigree(pedigrees[i]);
    }

    return ret;
}

bool LinkageProgram::run_pedigree(Pedigree& p) {
    unsigned iterations = 100000; //800 * p.num_members() * p.num_markers() * 10 * 2;
    SimwalkDescentGraph* opt;

/*
    opt = new SimwalkDescentGraph(&p, &map);
    opt->random_descentgraph();
    opt->likelihood();
*/
// <--- current
/*
    SimulatedAnnealing sa(&p, &map);
    opt = sa.optimise(iterations);
    
    Peeler peeler(&p, &map);
    peeler.peel(opt);
    
    delete opt;
    return true;
*/
/*    
    opt = new SimwalkDescentGraph(&p, &map);
    
    opt->random_descentgraph();
    opt->likelihood();
    printf("random descent graph prob = %f\n", opt->get_prob());
    opt->print();
    
    delete opt;
    exit(0);
*/    


    printf("processing pedigree %s\n", p.get_id().c_str());

    // run simulated annealing
    SimulatedAnnealing sa(&p, &map);
    opt = sa.optimise(iterations);
    printf("optimised prob = %f\n", opt->get_prob());
    
    // run markov chain
    MarkovChain mc(&p, &map);
    mc.run(opt, iterations);
    
    // print out results
    Peeler* peel = mc.get_peeler();
    peel->print();
    
        
    delete opt;
    
	return true; // XXX can this fail?

}

