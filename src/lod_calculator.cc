using namespace std;

#include <cmath>

#include "genetic_map.h"
#include "pedigree.h"
#include "lod_calculator.h"


LodCalculator::LodCalculator(Pedigree* p, GeneticMap* g) 
    : ped(p), map(g) {
    
    lod_scores = new double[map->num_markers() - 1];
    
    // XXX initialisation is a pain, really needs to be log(0) or -inf
    for(unsigned i = 0; i < (map->num_markers() - 1); ++i)
        lod_scores[i] = 0.0;
    
    random_prob = log10(1.0 / pow(2, float(ped->num_members() - ped->num_founders())));
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

void LodCalculator::add(unsigned locus, double prob) {
    lod_scores[locus] = log_add(
                            lod_scores[locus], 
                            log10(prob) - random_prob
                        );
}

double LodCalculator::get(unsigned locus) {
    return lod_scores[locus];
}

void LodCalculator::print() {
    printf("\nLOCUS\tLOD\n");
    for(unsigned i = 0; i < (map->num_markers() - 1); ++i) {
        printf("%d\t%f\n", i, lod_scores[i]);
    }
}

