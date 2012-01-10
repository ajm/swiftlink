#ifndef LKG_MARKOVCHAIN_H_
#define LKG_MARKOVCHAIN_H_

using namespace std;

class Pedigree;
class GeneticMap;
class Peeler;
class PeelSequenceGenerator;
class DescentGraph;


class MarkovChain {
    
    Pedigree* ped;
    GeneticMap* map;
    
    void initialise(DescentGraph& dg, PeelSequenceGenerator& psg);
    void parallel_initialise(DescentGraph& dg, PeelSequenceGenerator& psg);
    
 public :
    MarkovChain(Pedigree* ped, GeneticMap* map) :
        ped(ped), map(map) {}
    
    ~MarkovChain() {}
    
    MarkovChain(const MarkovChain& rhs) :
        ped(rhs.ped), map(map) {}
    
    MarkovChain& operator=(const MarkovChain& rhs) {
        
        if(this != &rhs) {
            ped = rhs.ped;
            map = rhs.map;
        }
        
        return *this;
    }
    
    double* run(unsigned iterations, double temperature);
};

#endif
