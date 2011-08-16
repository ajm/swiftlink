using namespace std;

#include <cstdio>
#include <cmath>

#include <iostream>
#include <sstream>
#include <string>

#include "genetic_map.h"


bool GeneticMap::sanity_check() {
    bool sane = ((map.size() - 1) == thetas.size()) and \
           ((map.size() - 1) == inversethetas.size());

    if(not sane) {
        fprintf(stderr, 
            "error in map data: number of markers = %d, number of thetas = %d\n", 
            int(map.size()), int(thetas.size()));
    }
    
    for(unsigned i = 0; i < map.size(); ++i)
        map[i].init_probs();
    
    for(unsigned i = 0; i < thetas.size(); ++i)
        half_thetas.push_back(haldane((get_marker(i+1).get_g_distance() - get_marker(i).get_g_distance()) / 2.0));
    
    for(unsigned i = 0; i < half_thetas.size(); ++i)
        half_inversethetas.push_back(1.0 - half_thetas[i]);
    
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

double GeneticMap::get_theta_halfway(unsigned int i) {
    return half_thetas[i];
}

double GeneticMap::get_inversetheta_halfway(unsigned int i) {
    return half_inversethetas[i];
}

double GeneticMap::get_theta(unsigned int i) {
    return thetas[i];
}

double GeneticMap::get_inversetheta(unsigned int i) {
    return inversethetas[i];
}

double GeneticMap::get_theta_log(unsigned int i) {
    return log(thetas[i]);
}

double GeneticMap::get_inversetheta_log(unsigned int i) {
    return log1p(-thetas[i]);
}

