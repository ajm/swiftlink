#include <cstdio>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <sstream>
#include <string>

#include "genetic_map.h"

using namespace std;


bool GeneticMap::sanity_check() {
    double tmp;

    if(map.size() != (thetas.size() + 1)) {
        bool dat_warning = false;

        fprintf(stderr, "Error: number of markers = %d, number of thetas = %d\n", int(map.size()), int(thetas.size()));

        for(int i = 0; i < int(map.size()); ++i) {
            if(not map[i].is_maf_set()) {
                fprintf(stderr, "Error: marker %d (%s) does not have a minor allele frequency\n", i+1, map[i].get_name().c_str());
                dat_warning = true;
            }
        }

        if(dat_warning) {
            fprintf(stderr, "Error: please check that DAT and MAP files are consistent...\n");
        }

        return false;
    }
    
    for(unsigned i = 0; i < map.size(); ++i)
        map[i].init_probs();
    
    // XXX dump map values from file, calculate directly from the map
    //     instead
    fprintf(stderr, "WARNING: recalculating theta values using Haldane map function\n");
    thetas.clear();
    inversethetas.clear();
    for(unsigned i = 1; i < map.size(); ++i) {
        tmp = haldane(map[i].get_g_distance() - map[i-1].get_g_distance());
        add_theta(tmp);

        if(tmp == 0.0) {
            fprintf(stderr, 
                    "Error: recombination fraction between %s and %s (markers %d and %d) is zero! (sampling will not work properly)\nExiting...\n", 
                    map[i-1].get_name().c_str(), map[i].get_name().c_str(), i-1, i);
            exit(EXIT_FAILURE);
        }
    }
    // XXX

    for(unsigned i = 0; i < thetas.size(); ++i) {
        tmp = haldane(inverse_haldane(thetas[i]) / double(partial_theta_count + 1));
        
        partial_thetas.push_back(tmp);
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

// coldest temperature is 1.0
// hottest is 0.0 (random)
void GeneticMap::set_temperature(double t) {
    if(temperature != 1.0) {
        fprintf(stderr, "error: temperature cannot be set twice in GeneticMap objects\n");
        abort();
    }
    
    temperature = t;
    
    // only theta and minor allele frequencies play 
    // a part in sampling, the partial theta stuff is only
    // used when calculating location scores
    for(unsigned i = 0; i < thetas.size(); ++i) {
        thetas[i] =         (temperature * thetas[i])         + ((1 - temperature) * 0.5);
        inversethetas[i] =  1.0 - thetas[i] ; //(temperature * inversethetas[i])  + ((1 - temperature) * 0.5);
    }

    for(unsigned i = 0; i < map.size(); ++i) {
        Snp& tmp = map[i];
        double minor_freq = tmp.minor();

        tmp.set_minor_freq((temperature * minor_freq) + ((1 - temperature) * 0.5));
    }
}

double GeneticMap::get_genetic_position(unsigned int index, unsigned int offset) const {
    return map[index].get_g_distance() + (offset ? inverse_haldane(partial_thetas[index] * offset) : 0);
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

