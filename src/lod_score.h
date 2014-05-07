#ifndef LKG_LODSCORE_H_
#define LKG_LODSCORE_H_

using namespace std;

#include <cstdio>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>

#include "genetic_map.h"
#include "logarithms.h"


class LODscores {
    
    GeneticMap* map;
    unsigned int num_scores_per_marker;
    unsigned int num_scores;
    unsigned int count;
    double trait_prob;
    vector<double> scores;
    vector<bool> initialised;
    
  public:
    LODscores(GeneticMap* map) : 
        map(map),
        num_scores_per_marker(map->get_lodscore_count()),
        num_scores((num_scores_per_marker * (map->num_markers() - 1))),
        count(0),
        trait_prob(0.0),
        scores(num_scores),
        initialised(num_scores) {}
        
    ~LODscores() {}
    
    LODscores(const LODscores& rhs) :
        map(rhs.map),
        num_scores_per_marker(rhs.num_scores_per_marker),
        num_scores(rhs.num_scores),
        count(rhs.count),
        trait_prob(rhs.trait_prob),
        scores(rhs.scores),
        initialised(rhs.initialised) {}
    
    LODscores& operator=(const LODscores& rhs) {
        
        if(&rhs != this) {
            map = rhs.map;
            num_scores_per_marker = rhs.num_scores_per_marker;
            num_scores = rhs.num_scores;
            count = rhs.count;
            trait_prob = rhs.trait_prob;
            scores = rhs.scores;
            initialised = rhs.initialised;
        }
        
        return *this;
    }
    
    void set_trait_prob(double prob) {
        trait_prob = prob;
    }
    
    unsigned int num_lodscores() const {
        return num_scores_per_marker * (map->num_markers() - 1);
    }
    
    unsigned int get_lodscores_per_marker() const {
        return num_scores_per_marker;
    }
    
    void add(unsigned int locus, unsigned int offset, double prob) {
        unsigned int index = (locus * num_scores_per_marker) + offset;
        scores[index] = initialised[index] ? log_sum(prob, scores[index]) : prob, initialised[index] = true;
        
        if((locus == 0) and (offset == 0))
            ++count;
    }

    double get_raw(unsigned int index) {
        return scores[index];
    }
    
    double get(unsigned int locus, unsigned int offset) const {
        return (scores[(locus * num_scores_per_marker) + offset] - log(count) - trait_prob) / log(10.0);
    }
    
    double get_genetic_position(unsigned int locus, unsigned int offset) {
        return map->get_genetic_position(locus, offset);
    }

    unsigned int get_count() {
        return count;
    }

    void merge_results(LODscores* tmp) {
        for(unsigned i = 0; i < map->num_markers() - 1; ++i) {
            for(unsigned j = 0; j < num_scores_per_marker; ++j) {
                unsigned int index = (i * num_scores_per_marker) + j;
                scores[index] = log_sum(scores[index], tmp->get_raw(index));
            }
        }
        count += tmp->get_count();
    }
    
    void set_count(unsigned int c) { count = c; }
    void set(unsigned int index, double prob) {
        //unsigned int index = (locus * num_scores_per_marker) + offset;
        
        scores[index] = prob;
        initialised[index] = true;
    }

    string debug_string() {
        stringstream ss;

        ss << "LodScore: ";

        for(unsigned int i = 0; i < scores.size(); ++i)
            ss << setprecision(2) << scores[i] << " ";

        return ss.str();
    }
};

#endif

