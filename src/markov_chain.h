#ifndef LKG_MARKOVCHAIN_H_
#define LKG_MARKOVCHAIN_H_

using namespace std;

#include "types.h"
#include "genetic_map.h"

class Pedigree;
class PeelSequenceGenerator;
class DescentGraph;
class LODscores;


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
    
        //map->set_temperature(0.0);    
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
    
    LODscores* run(DescentGraph& dg);
    bool noninterferring(vector<int>& x, int val);
};

#endif

