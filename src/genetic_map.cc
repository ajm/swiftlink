using namespace std;

#include <cstdio>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <sstream>
#include <string>

#include "genetic_map.h"


bool GeneticMap::sanity_check() {
    double tmp;

    if(map.size() != (thetas.size() + 1)) {
        fprintf(stderr, "error: number of markers = %d, number of thetas = %d\n", int(map.size()), int(thetas.size()));
        return false;
    }
    
    for(unsigned i = 0; i < map.size(); ++i)
        map[i].init_probs();
    
    for(unsigned i = 0; i < thetas.size(); ++i) {
        tmp = haldane(inverse_haldane(thetas[i]) / double(partial_theta_count + 1));
        
        partial_thetas.push_back(tmp);
    }
    
    exit(0);
    
    return true;
}

string GeneticMap::debug_string() {
    stringstream ss;
    
    ss.precision(DEBUG_FP_PRECISION);
    
    ss << "GeneticMap: " << map.size() << " loci" << endl;
    
    for(unsigned i = 0; i < map.size(); ++i) {
        ss << scientific << map[i].debug_string() << endl;
    }
    
    ss << "  thetas:" << endl;
    for(unsigned i = 0; i < thetas.size(); ++i) {
        ss << "\t" << scientific << exp(thetas[i]) << endl;
    }
    
    return ss.str();
}

double GeneticMap::haldane(double m) const {
    return 0.5 * (1.0 - exp(-2.0 * m));
}

double GeneticMap::inverse_haldane(double r) const {
    return -0.5 * log(1 - (2 * r));
}

// XXX
// this was never used, the idea was to do parallel tempering
// but in the end it did not seem necessary
void GeneticMap::set_temperature(double t) {
    if(temperature != 0.0) {
        fprintf(stderr, "error: temperature cannot be set twice in GeneticMap objects\n");
        abort();
    }
    
    temperature = t;
    
    for(unsigned i = 0; i < thetas.size(); ++i) {
        thetas[i] =         ((1 - temperature) * thetas[i])         + (temperature * 0.5);
        inversethetas[i] =  ((1 - temperature) * inversethetas[i])  + (temperature * 0.5);
    }
}

double GeneticMap::get_genetic_position(unsigned int index, unsigned int offset) const {
    return map[index].get_g_distance() + inverse_haldane(partial_thetas[index] * offset);
}

// XXX offset is 1 -- partial_theta_count
double GeneticMap::get_theta_partial(unsigned int index, unsigned int offset) const {
    return partial_thetas[index] * offset;
}

double GeneticMap::get_inversetheta_partial(unsigned int index, unsigned int offset) const {
    return 1.0 - (partial_thetas[index] * offset);
}

double GeneticMap::get_theta_log(unsigned int i) const {
    return log(thetas[i]);
}

double GeneticMap::get_inversetheta_log(unsigned int i) const {
    return log(inversethetas[i]);
}

