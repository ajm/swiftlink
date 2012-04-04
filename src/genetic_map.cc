using namespace std;

#include <cstdio>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <sstream>
#include <string>

#include "genetic_map.h"


bool GeneticMap::sanity_check() {
    if(map.size() != (thetas.size() + 1)) {
        fprintf(stderr, "error: number of markers = %d, number of thetas = %d\n", int(map.size()), int(thetas.size()));
        return false;
    }
    
    for(unsigned i = 0; i < map.size(); ++i)
        map[i].init_probs();
    
    for(unsigned i = 0; i < thetas.size(); ++i) {
        //double tmp = haldane((get_marker(i+1).get_g_distance() - get_marker(i).get_g_distance()) / 2.0);
        
        // XXX I am not sure if I can trust the genetic distances in the map
        // file, this way at least I am consistent irrespective of what the user
        // states the map is...
        
        double tmp = haldane(inverse_haldane(thetas[i]) / 2.0);
        
        half_thetas.push_back(tmp);
        half_inversethetas.push_back(1.0 - tmp);
    }
    
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

double GeneticMap::get_theta_halfway(unsigned int i) const {
    return half_thetas[i];
}

double GeneticMap::get_inversetheta_halfway(unsigned int i) const {
    return half_inversethetas[i];
}
/*
double GeneticMap::get_theta(unsigned int i) const {
    return thetas[i];
}

double GeneticMap::get_inversetheta(unsigned int i) const {
    return inversethetas[i];
}
*/
double GeneticMap::get_theta_log(unsigned int i) const {
    return log(thetas[i]);
}

double GeneticMap::get_inversetheta_log(unsigned int i) const {
    return log(inversethetas[i]);
}

