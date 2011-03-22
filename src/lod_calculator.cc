using namespace std;

#include <cmath>
#include <limits>

#include "genetic_map.h"
#include "pedigree.h"
#include "lod_calculator.h"


LodCalculator::LodCalculator(Pedigree* p, GeneticMap* g) 
    : ped(p), map(g), initialised(g->num_markers() - 1, false), count(0) {
    
    lod_scores = new double[map->num_markers() - 1];
    random_prob = log10(1.0 / pow(2, double(ped->num_members() - ped->num_founders())));
}

LodCalculator::~LodCalculator() {
    delete[] lod_scores;
}

inline double LodCalculator::exp10(double x) {
    return pow(10.0, x);
}

double LodCalculator::log_add(double a, double b) {
    // log(exp(a - b) + exp(b - b)) + b
    return log10(exp10(a - b) + 1) + b;
}

void LodCalculator::add(unsigned locus, double prob, double trans) {
/*    
    if(locus == 0)
        printf("locus = %d, prob = %e\n", locus, prob);
*/
    lod_scores[locus] = initialised[locus] ? 
        log_add(lod_scores[locus], log10(prob) - trans /*random_prob*/) : 
        log10(prob) - trans /*random_prob*/, initialised[locus] = true;
    if(locus == 0)
        count++;
}

double LodCalculator::get(unsigned locus) {
    return lod_scores[locus];
}

void LodCalculator::print() {
    double tot = log10(count);
    printf("count = %d (tot = %e)\n", count, tot);
    printf("\nLOCUS\tLOD\n");
    for(unsigned i = 0; i < (map->num_markers() - 1); ++i) {
        printf("%d\t%f\n", 
            i, lod_scores[i] == -numeric_limits<double>::infinity() ? 
                    lod_scores[i] : 
                    (lod_scores[i] - tot) - random_prob);
    }
}

