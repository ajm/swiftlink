#ifndef LKG_SEQUENTIALIMPUTATION_H_
#define LKG_SEQUENTIALIMPUTATION_H_

using namespace std;


class Pedigree;
class GeneticMap;
class PeelSequenceGenerator;
class DescentGraph;


class SequentialImputation {
    
    Pedigree* ped;
    GeneticMap* map;
    PeelSequenceGenerator* psg;
    
 public :
    SequentialImputation(Pedigree* ped, GeneticMap* map, PeelSequenceGenerator* psg) :
        ped(ped), 
        map(map),
        psg(psg) {}
    
    ~SequentialImputation() {}
    
    SequentialImputation(const SequentialImputation& rhs) :
        ped(rhs.ped), 
        map(rhs.map),
        psg(rhs.psg) {}
    
    SequentialImputation& operator=(const SequentialImputation& rhs) {
        if(this != &rhs) {
            ped = rhs.ped;
            map = rhs.map;
            psg = rhs.psg;
        }
        return *this;
    }
    
    void run(DescentGraph& dg, int iterations);
    void parallel_run(DescentGraph& dg, int iterations);
};

#endif

