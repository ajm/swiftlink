#ifndef LKG_MC3_
#define LKG_MC3_

using namespace std;

#include <vector>

#include "types.h"

class Pedigree;
class GeneticMap;
class PeelSequenceGenerator;
class MarkovChain;
class LODscores;

class Mc3 {
    
    Pedigree* ped;
    GeneticMap* map;
    struct mcmc_options options;

    PeelSequenceGenerator* psg;
    LODscores* lod;
    vector<Peeler*> peelers;
    vector<MarkovChain*> chains;

    void _init();
    void _kill();

  public:
    Mc3(Pedigree* ped, GeneticMap* map, PeelSequenceGenerator* psg, struct mcmc_options opt) :
        ped(ped),
        map(map),
        options(opt),
        psg(psg),
        lod(0),
        peelers(), 
        chains() {
    
        _init();
    }

    ~Mc3() {
        _kill();
    }

    Mc3(const Mc3& rhs) :
        ped(rhs.ped),
        map(rhs.map),
        options(rhs.options),
        psg(rhs.psg),
        lod(rhs.lod),
        peelers(rhs.peelers),
        chains(rhs.chains) {}

    Mc3& operator=(const Mc3& rhs) {
        if(this != &rhs) {
            ped = rhs.ped;
            map = rhs.map;
            options = rhs.options;
            psg = rhs.psg;
            lod = rhs.lod;
            peelers = rhs.peelers;
            chains = rhs.chains;
        }
        return *this;
    }

    LODscores* run(DescentGraph& dg);
};

#endif

