#ifndef LKG_PEELER_H_
#define LKG_PEELER_H_

using namespace std;

#include <vector>

#include "peeling.h"
#include "rfunction.h"


class Pedigree;
class GeneticMap;
class SimwalkDescentGraph;

class Peeler {
    
    Pedigree* ped;
    GeneticMap* map;
    vector<Rfunction> rfunctions;
        
 public :
    Peeler(Pedigree* p, GeneticMap* g);

    bool peel(SimwalkDescentGraph* sdg);
};

#endif

