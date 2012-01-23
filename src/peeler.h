#ifndef LKG_PEELER_H_
#define LKG_PEELER_H_

using namespace std;

#include <vector>

#include "peeling.h"
#include "trait_rfunction.h"
#include "peel_sequence_generator.h"


class Pedigree;
class GeneticMap;
class DescentGraph;

class Peeler {
    
    Pedigree* ped;
    GeneticMap* map;
    vector<TraitRfunction> rfunctions;
    double trait_prob;
    double lod_score;
    bool initialised;
    unsigned int locus;
    unsigned int count;
    
    double calc_trait_prob();
    
 public :
    Peeler(Pedigree* p, GeneticMap* g, PeelSequenceGenerator* psg, unsigned int locus);
    Peeler(const Peeler& rhs);
    ~Peeler();
    
    Peeler& operator=(const Peeler& rhs);
    
    double get_trait_prob();
    void process(DescentGraph* dg);
    double get();
};

#endif
