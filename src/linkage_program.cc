using namespace std;

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "descent_graph.h"
#include "simulated_annealing.h"
#include "linkage_program.h"
#include "linkage_writer.h"
#include "disease_model.h"
#include "markov_chain.h"
#include "genetic_map.h"
#include "pedigree.h"
#include "peeler.h"
#include "locus_sampler.h"


bool LinkageProgram::run() {
    bool ret = true;
    
    if(verbose) {
        fprintf(stderr, "%s\n", dm.debug_string().c_str());
        fprintf(stderr, "%s\n", map.debug_string().c_str());
    }

    // TODO XXX need to know how to do this properly, 
    // look up better random numbers for simulations etc
    srandom(time(NULL));

    for(unsigned int i = 0; i < pedigrees.size(); ++i) {
        if(verbose) {
            fprintf(stderr, "%s\n", pedigrees[i].debug_string().c_str());
        }
        
        ret &= run_pedigree(pedigrees[i]);
    }

    return ret;
}

bool LinkageProgram::run_pedigree(Pedigree& p) {
    // TODO XXX this should probably go somewhere else/be user-defined
//    unsigned iterations = 1600 * p.num_members() * p.num_markers() * 20 * 2;
//    DescentGraph* opt;
    Peeler* peel;    
    
    if(verbose) {
        fprintf(stderr, "processing pedigree %s\n", p.get_id().c_str());
    }
/*    
    // run simulated annealing
    SimulatedAnnealing sa(&p, &map);
    opt = sa.optimise(iterations);
    
    // run markov chain
    MarkovChain mc(&p, &map);
    peel = mc.run(opt, iterations);
*/

    LocusSampler lsampler(&p, &map);
    //peel = lsampler.run(10000);
    peel = lsampler.temper(100000, 10);
    
    // write out results
    LinkageWriter lw(&map, peel, "linkage.txt", verbose);
    lw.write();

    // TODO XXX I should not write out immediately, but store the results
    // combine them and then write out everything in a table
    
//    delete opt;
    
	return true;
}

