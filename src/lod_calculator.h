#ifndef LKG_LODCALCULATOR_H_
#define LKG_LODCALCULATOR_H_

using namespace std;


class Pedigree;
class GeneticMap;

class LodCalculator {
    
    Pedigree* ped;
    GeneticMap* map;
    double* lod_scores;
    double random_prob;
    
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

