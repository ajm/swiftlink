#ifndef LKG_MARKOVCHAIN_H_
#define LKG_MARKOVCHAIN_H_

using namespace std;

#include "types.h"
#include "genetic_map.h"

class Pedigree;
class PeelSequenceGenerator;
class DescentGraph;


class MarkovChain {
    
    Pedigree* ped;
    GeneticMap* map;
    PeelSequenceGenerator* psg;
    struct mcmc_options options;
    
 public :
    MarkovChain(Pedigree* ped, GeneticMap* map, PeelSequenceGenerator* psg, struct mcmc_options options) :
        ped(ped), 
        map(map), 
        psg(psg),
        options(options) {
    
        map->set_temperature(options.temperature);    
    }
    
    ~MarkovChain() {}
    
    MarkovChain(const MarkovChain& rhs) :
        ped(rhs.ped), 
        map(rhs.map),
        psg(rhs.psg), 
        options(rhs.options) {}
    
    MarkovChain& operator=(const MarkovChain& rhs) {
        if(this != &rhs) {
            ped = rhs.ped;
            map = rhs.map;
            psg = rhs.psg;
            options = rhs.options;
        }
        return *this;
    }
    
    double* run(DescentGraph& dg);
};

#endif

