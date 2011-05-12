using namespace std;

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>

#include "haplotype_program.h"
#include "genetic_map.h"
#include "disease_model.h"
#include "descent_graph.h"
#include "simwalk_descent_graph.h"
#include "simulated_annealing.h"
#include "haplotype_writer.h"
#include "pedigree.h"


bool HaplotypeProgram::run() {
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
        
        ret &= run_pedigree_sa(pedigrees[i]);
    }

    return ret;
}

bool HaplotypeProgram::run_pedigree_sa(Pedigree& p) {
    //unsigned iterations = 800 * p.num_members() * p.num_markers() * 10 * 2;
    unsigned iterations = 100000; // for testing
    SimwalkDescentGraph* opt;
        
    if(verbose) {
        fprintf(stderr, "processing pedigree %s\n", p.get_id().c_str());
    }

    // run simulated annealing
    // XXX TODO this is using the wrong objective function,
    // it is optimising the descent graph likelihood, not the 
    // descent state likelihood 
    SimulatedAnnealing sa(&p, &map);
    opt = sa.optimise(iterations);
    
    stringstream ss;
    ss << "haplotype_" << p.get_id() << ".txt";
    
    HaplotypeWriter hw(static_cast<DescentGraph*>(opt), &p, &map, ss.str().c_str(), verbose);
    hw.write();
    
    delete opt;

    return true;
}

