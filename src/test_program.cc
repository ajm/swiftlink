using namespace std;

#include <cstdio>

#include "test_program.h"
#include "descent_graph.h"
#include "disease_model.h"
#include "genetic_map.h"
#include "pedigree.h"
#include "peel_sequence_generator.h"
//#include "gpu_wrapper.h"

bool TestProgram::run() {

    Pedigree* ped = &(pedigrees[0]);

    PeelSequenceGenerator psg(ped);
    psg.build_peel_order();

    //GPUWrapper gpu(ped, &map, psg);
    
    fprintf(stderr, "test program not used at the moment...\n");
    abort();
    
    return true;
}

