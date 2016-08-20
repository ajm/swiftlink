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
class LODscores;

class Peeler {
    
    Pedigree* ped;
    GeneticMap* map;
    LODscores* lod;
    vector<TraitRfunction> rfunctions;
    unsigned int locus;
    bool sex_linked;
    
 public :
    Peeler(Pedigree* p, GeneticMap* g, PeelSequenceGenerator* psg, LODscores* lod, bool sex_linked);
    Peeler(const Peeler& rhs);
    ~Peeler();
    
    Peeler& operator=(const Peeler& rhs);
    
    double calc_trait_prob();
    double get_trait_prob(); // XXX to ease transition, but delete later...
    
    void set_locus(unsigned int l) {
        locus = l;
        
        for(unsigned i = 0; i < rfunctions.size(); ++i) {
            rfunctions[i].set_locus_minimal(locus);
        }
    }
    
    void process(DescentGraph* dg);
};

#endif

