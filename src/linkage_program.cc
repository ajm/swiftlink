using namespace std;

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "descent_graph.h"
#include "linkage_program.h"
#include "linkage_writer.h"
#include "disease_model.h"
#include "markov_chain.h"
#include "genetic_map.h"
#include "pedigree.h"
#include "peeler.h"


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
    
    if(verbose) {
        fprintf(stderr, "processing pedigree %s\n", p.get_id().c_str());
    }

    MarkovChain chain(&p, &map);
    //Peeler* peel = chain.run(1000000, 0.0);
    Peeler* peel = chain.run(100, 0.0);
    
    // write out results
    LinkageWriter lw(&map, peel, "linkage.txt", verbose);
    lw.write();

    // TODO XXX I should not write out immediately, but store the results
    // combine them and then write out everything in a table

    delete peel;
    
    return true;
}
