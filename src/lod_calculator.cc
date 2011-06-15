using namespace std;

#include <cmath>
#include <limits>

#include "genetic_map.h"
#include "pedigree.h"
#include "lod_calculator.h"


LodCalculator::LodCalculator(Pedigree* p, GeneticMap* g) : 
    ped(p), 
    map(g), 
    lod_scores(NULL),
    initialised(map->num_markers() - 1, false), 
    count(0), 
    trait_prob(0.0) {
    
    lod_scores = new double[map->num_markers() - 1];
}

LodCalculator::LodCalculator(const LodCalculator& l) :
    ped(l.ped),
    map(l.map),
    lod_scores(NULL),
    initialised(l.initialised),
    count(l.count),
    trait_prob(l.trait_prob) {

    lod_scores = new double[map->num_markers() - 1];
    copy(l.lod_scores, 
         l.lod_scores + (map->num_markers() - 1), 
         lod_scores);
}

LodCalculator::~LodCalculator() {
    delete[] lod_scores;
}

LodCalculator& LodCalculator::operator=(const LodCalculator& rhs) {

    if(&rhs != this) {
        if(map->num_markers() != rhs.map->num_markers()) {
            delete[] lod_scores;
            lod_scores = new double[map->num_markers() - 1];
        }
    
        ped = rhs.ped;
        map = rhs.map;
        initialised = rhs.initialised;
        count = rhs.count;
        trait_prob = rhs.trait_prob;
        
        copy(rhs.lod_scores, 
             rhs.lod_scores + (map->num_markers() - 1), 
             lod_scores);
    }
    
    return *this;
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

