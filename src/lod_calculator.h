#ifndef LKG_LODCALCULATOR_H_
#define LKG_LODCALCULATOR_H_

using namespace std;

#include <vector>


class Pedigree;
class GeneticMap;

class LodCalculator {
    
    Pedigree* ped;
    GeneticMap* map;
    double* lod_scores;
    vector<bool> initialised; // XXX <-- I don't like this, but log(0.0) is -inf :-S
    unsigned int count;
    double trait_prob;
        
 public :
    LodCalculator(Pedigree* p, GeneticMap* g);
    LodCalculator(const LodCalculator& l);
    ~LodCalculator();
    LodCalculator& operator=(const LodCalculator& rhs);
    
    void add(unsigned locus, double prob);
    double get(unsigned locus);
    
    void set_trait_prob(double p);
};

#endif

