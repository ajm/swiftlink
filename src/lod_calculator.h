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
    double random_prob;
    unsigned count;
    
    double log_add(double a, double b);
    double exp10(double x);
    
 public :
    LodCalculator(Pedigree* p, GeneticMap* g);
    ~LodCalculator();
    
    void add(unsigned locus, double prob);
    double get(unsigned locus);
    void print();
};

#endif

