using namespace std;

#include <cmath>

#include "markov_chain.h"
#include "peel_sequence_generator.h"
#include "peeler.h"
#include "locus_sampler2.h"
#include "meiosis_sampler.h"
#include "pedigree.h"
#include "genetic_map.h"
#include "descent_graph.h"
#include "linkage_writer.h"
#include "progress.h"


Peeler* MarkovChain::run(unsigned iterations, double temperature) {
    unsigned burnin = iterations * 0.01;
    
    map->set_temperature(temperature);

    // create a descent graph
    DescentGraph dg(ped, map);
    dg.random_descentgraph();
    
    /*
    // start at an optimal configuration for east
    dg.set(14, 0, PATERNAL, 1);
    dg.set(14, 1, PATERNAL, 1);
    dg.set(14, 2, PATERNAL, 1);
    dg.set(14, 3, PATERNAL, 1);
    dg.set(14, 4, PATERNAL, 1);
    dg.set(14, 5, PATERNAL, 1);
    dg.set(14, 6, PATERNAL, 1);
    dg.set(14, 7, PATERNAL, 1);
    dg.set(14, 8, PATERNAL, 1);
    dg.set(14, 9, PATERNAL, 1);
    */
    
    // build peeling sequence for L-sampler and Peeler
    PeelSequenceGenerator psg(ped);
    psg.build_peel_order();
    
    // create an object to perform LOD scoring
    // allocate on the heap so we can return it
    Peeler* peel = new Peeler(ped, map, psg);
    
    // create samplers
    LocusSampler lsampler(ped, map, psg);
    MeiosisSampler msampler(ped, map);
    
    // some vars to keep track of what is being sampled
    // just cycle through markers and people
    unsigned locus = 0;
    unsigned person = ped->num_founders(); // this will be the index of the first non-founder

    //Progress p("MCMC: ", iterations);
    
    for(unsigned i = 0; i < iterations; ++i) {
        //if((i % 100) == 0)
        //    printf("%d\n", i);
        
        //lsampler.step(dg, locus);
        //locus = (locus + 1) % map->num_markers();
        
        if((random() / static_cast<double>(RAND_MAX)) < 0.5) {
            //printf("L");
            lsampler.step(dg, locus);
            locus = (locus + 1) % map->num_markers();
        }
        else {
            //printf("M");
            msampler.step(dg, person);
            person = (person + 1) % ped->num_members();
            if(person == 0) {
                person = ped->num_founders();
            }
        }
        
        //printf("X %f %d\n", dg._recombination_prob(), dg.num_recombinations());
        
        //p.increment();
        
        if(i < burnin) {
            continue;
        }
        
        if((i % 10) == 0)
            peel->process(dg);
        
        
        if((i % 10) == 0) {
            LinkageWriter lw(map, peel, "linkage.txt", true);
            lw.write();
        }
        
    }
    
    //p.finish();
    
    return peel;
}
