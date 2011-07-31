using namespace std;

#include <cstdio>
#include <cmath>

#include <iostream>
#include <sstream>
#include <string>

#include "genetic_map.h"


bool GeneticMap::sanity_check() {
    bool sane = ((map.size() - 1) == thetas.size()) and \
           ((map.size() - 1) == inverse_thetas.size());

    if(not sane) {
        fprintf(stderr, 
            "error in map data: number of markers = %d, number of thetas = %d\n", 
            int(map.size()), int(thetas.size()));
    }
    
    return sane;
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

double GeneticMap::haldane(double m) {
    return 0.5 * (1.0 - exp(-2.0 * m));
}

double GeneticMap::get_theta_halfway(unsigned int i) {
    return haldane((get_marker(i+1).get_g_distance() - get_marker(i).get_g_distance()) / 2.0);
}

double GeneticMap::get_theta(unsigned int i) {
    return log(((1 - temperature) * exp(thetas[i])) + (temperature * 0.5));
}

double GeneticMap::get_inverse_theta(unsigned int i) {
    return log(((1 - temperature) * exp(inverse_thetas[i])) + (temperature * 0.5));
}

