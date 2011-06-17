#ifndef LKG_PEELER_H_
#define LKG_PEELER_H_

using namespace std;

#include <vector>

#include "peeling.h"
#include "rfunction.h"
#include "lod_calculator.h"


class Pedigree;
class GeneticMap;
class DescentGraph;

class Peeler {
    
    Pedigree* ped;
    GeneticMap* map;
    vector<Rfunction*> rfunctions;
    LodCalculator lod;
    double trait_prob;
    
    
    double calc_trait_prob();
    void copy_rfunctions(const Peeler& rhs);
    
 public :
    Peeler(Pedigree* p, GeneticMap* g);
    Peeler(const Peeler& rhs);
    ~Peeler();
    
    Peeler& operator=(const Peeler& rhs);
    
    double peel(DescentGraph* dg, unsigned locus);
    double get_trait_prob();
    bool process(DescentGraph& dg);
    double get(unsigned locus);
};

#endif

