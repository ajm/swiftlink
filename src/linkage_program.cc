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
    vector<Peeler*> peelers;
    Peeler* tmp;
    
    if(verbose) {
        //fprintf(stderr, "%s\n", dm.debug_string().c_str());
        //fprintf(stderr, "%s\n", map.debug_string().c_str());
    }

    // TODO XXX need to know how to do this properly, 
    // look up better random numbers for simulations etc
    srandom(time(NULL));

    for(unsigned int i = 0; i < pedigrees.size(); ++i) {
        if(verbose) {
            fprintf(stderr, "%s\n", pedigrees[i].debug_string().c_str());
        }
        
        // it cannot actually be NULL, the program will call
        // abort() at the slightest hint of a problem
        if((tmp = run_pedigree(pedigrees[i])) == NULL) {
            fprintf(stderr, "error: pedigree '%s' failed\n", pedigrees[i].get_id().c_str());
            free_peelers(peelers);
            return false;
        }
        
        peelers.push_back(tmp);
    }
    
    LinkageWriter lw(&map, output_filename, verbose);
    if(not lw.write(peelers)) {
        fprintf(stderr, "error: could not write output file '%s'\n", output_filename.c_str());
        return false;
    }
    
    free_peelers(peelers);

    return true;
}

Peeler* LinkageProgram::run_pedigree(Pedigree& p) {
    
    if(verbose) {
        fprintf(stderr, "processing pedigree %s\n", p.get_id().c_str());
    }

    MarkovChain chain(&p, &map);
    
    return chain.run(iterations, 0.0);
}

void LinkageProgram::free_peelers(vector<Peeler*>& p) {
    for(unsigned i = 0; i < p.size(); ++i)
        delete p[i];
}

