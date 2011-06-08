using namespace std;

#include <cmath>
#include <limits>

#include "genetic_map.h"
#include "pedigree.h"
#include "lod_calculator.h"


LodCalculator::LodCalculator(Pedigree& p, GeneticMap& g) 
    : ped(p), 
      map(g), 
      initialised(g.num_markers() - 1, false), 
      count(0), 
      trait_prob(0.0) {
    
    lod_scores = new double[map.num_markers() - 1];
}

LodCalculator::~LodCalculator() {
    delete[] lod_scores;
}

void LodCalculator::set_trait_prob(double p) {
    trait_prob = p;
}

inline double LodCalculator::exp10(double x) {
    return pow(10.0, x);
}

double LodCalculator::log_add(double a, double b) {
    // log(exp(a - b) + exp(b - b)) + b
    return log10(exp10(a - b) + 1) + b;
}

void LodCalculator::add(unsigned locus, double prob) {

    lod_scores[locus] = initialised[locus] ? 
        log_add(lod_scores[locus], prob) : 
        prob, initialised[locus] = true;

    if(locus == 0)
        count++;
}

double LodCalculator::get(unsigned locus) {
    return (lod_scores[locus] - log10(count) - trait_prob) / log(10.0);
}

