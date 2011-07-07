#ifndef LKG_PARALLELTEMPER_H_
#define LKG_PARALLELTEMPER_H_

using namespace std;

#include <vector>

#include "locus_sampler.h"


class Peeler;
class Pedigree;
class GeneticMap;

class ParallelTempering {

    vector<LocusSampler*> chains;
    vector<double> temperatures;
    Pedigree* ped;
    GeneticMap* map;
    Peeler peel;
    
    void init_chains(unsigned num_chains);
    void copy_chains(const ParallelTempering& rhs);
    void kill_chains();
    bool exchange_replicas(LocusSampler* ls1, LocusSampler* ls2, double temp1, double temp2);
    double get_random();
    
 public :
    ParallelTempering(Pedigree* ped, GeneticMap* map, unsigned num_chains) :
        chains(),
        temperatures(),
        ped(ped),
        map(map),
        peel(ped, map) {
        
        init_chains(num_chains);
    }
    
    ParallelTempering(const ParallelTempering& rhs) :
        chains(),
        temperatures(rhs.temperatures),
        ped(rhs.ped),
        map(rhs.map),
        peel(rhs.peel) {
    
        copy_chains(rhs);
    }
    
    ParallelTempering& operator=(const ParallelTempering& rhs) {
        
        if(&rhs != this) {
            temperatures = rhs.temperatures;
            ped = rhs.ped;
            map = rhs.map;
            peel = rhs.peel;
            
            kill_chains();
            copy_chains(rhs);
        }
        
        return *this;
    }
    
    ~ParallelTempering() {
        kill_chains();
    }
    
    Peeler* run(unsigned iterations);
};

#endif

