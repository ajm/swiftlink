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
    bool sex_linked;
    
 public :
    SequentialImputation(Pedigree* ped, GeneticMap* map, PeelSequenceGenerator* psg, bool sex_linked) :
        ped(ped), 
        map(map),
        psg(psg),
        sex_linked(sex_linked) {}
    
    ~SequentialImputation() {}
    
    SequentialImputation(const SequentialImputation& rhs) :
        ped(rhs.ped), 
        map(rhs.map),
        psg(rhs.psg),
        sex_linked(rhs.sex_linked) {}
    
    SequentialImputation& operator=(const SequentialImputation& rhs) {
        if(this != &rhs) {
            ped = rhs.ped;
            map = rhs.map;
            psg = rhs.psg;
            sex_linked = rhs.sex_linked;
        }
        return *this;
    }
    
    void run(DescentGraph& dg, int iterations);
    void parallel_run(DescentGraph& dg, int iterations);
};

#endif

