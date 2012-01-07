#ifndef LKG_PEELER_H_
#define LKG_PEELER_H_

using namespace std;

#include <vector>

#include "peeling.h"
#include "trait_rfunction.h"
#include "lod_calculator.h"
#include "peel_sequence_generator.h"


class Pedigree;
class GeneticMap;
class DescentGraph;

class Peeler {
    
    Pedigree* ped;
    GeneticMap* map;
    vector<TraitRfunction> rfunctions;
    LodCalculator lod;
    double trait_prob;
    
    
    double calc_trait_prob();
    //void copy_rfunctions(const Peeler& rhs);
    //void kill_rfunctions();
    
 public :
    Peeler(Pedigree* p, GeneticMap* g, PeelSequenceGenerator& psg);
    Peeler(const Peeler& rhs);
    ~Peeler();
    
    Peeler& operator=(const Peeler& rhs);
    
    double peel(DescentGraph* dg, unsigned locus);
    double get_trait_prob();
    void process(DescentGraph& dg);
    double get(unsigned locus);
};

#endif

