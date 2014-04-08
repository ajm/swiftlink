#ifndef LKG_GENETICMAP_H_
#define LKG_GENETICMAP_H_

using namespace std;

#include <cstdio>
#include <cmath>

#include <vector>
#include <iostream>
#include <sstream>
#include <string>

#include "types.h"


// XXX create a super class 'marker'?
// or maybe I will just ignore more polymorphic markers
class Snp {

    string name;
    double genetic_distance;
    unsigned int physical_distance;
    double major_freq;
    double minor_freq;
    double prob[4];

 public :
    Snp(string name, double genetic, unsigned int physical=0) : 
        name(name), 
        genetic_distance(genetic), 
        physical_distance(physical),
        major_freq(1.0), 
        minor_freq(0.0) {}

    ~Snp() {}
    
    double major() const { return major_freq; }
    double minor() const { return minor_freq; }
    void set_minor_freq(double m) {
        minor_freq = m;
        major_freq = 1.0 - m;
    }
    
    string get_name() const { return name; }
    double get_g_distance() const { return genetic_distance; }
    unsigned int get_p_distance() const { return physical_distance; }
	
	string debug_string() {
	    stringstream ss;
	    
	    ss.precision(DEBUG_FP_PRECISION);
	    
	    ss << "\t" << name \
	       << "\t" << fixed << genetic_distance \
	       << "\t" << physical_distance \
	       << "\t(minor=" << fixed << minor_freq \
	       << ", major=" << fixed << major_freq << ")";
	    
	    return ss.str();
	}
	
	void init_probs() {
	    prob[TRAIT_UU] = major_freq * major_freq;
	    prob[TRAIT_AU] = prob[TRAIT_UA] = minor_freq * major_freq;
	    prob[TRAIT_AA] = minor_freq * minor_freq;
	    
	    double total = prob[TRAIT_UU] + prob[TRAIT_AU] + prob[TRAIT_UA] + prob[TRAIT_AA];
	    
	    for(int i = 0; i < 4; ++i)
	        prob[i] /= total;
	}
	
	double get_prob(enum phased_trait pt) const {
	    return prob[pt];
	}
};

class GeneticMap {

    vector<Snp> map;
    vector<double> thetas;
    vector<double> inversethetas;
    vector<double> partial_thetas;
    double temperature;
    unsigned int partial_theta_count; // must be greater than zero

    //double haldane(double m) const ;
    //double inverse_haldane(double m) const;
    
 public :
    GeneticMap(unsigned int partial_theta_count) : 
        map(), 
        thetas(), 
        inversethetas(),
        partial_thetas(),
        temperature(0.0),
        partial_theta_count(partial_theta_count) {}
    
    ~GeneticMap() {}
    
	Snp& operator[](int i) {
		return map[i];
	}
    
    const Snp& operator[](int i) const {
		return map[i];
	}

    double haldane(double m) const ;
    double inverse_haldane(double m) const;
    
    void add(Snp& s) {
        map.push_back(s);
    }
    
    void add_theta(double d) {
        //thetas.push_back(log(d));
        //inverse_thetas.push_back(log1p(-d));
        thetas.push_back(d);
        inversethetas.push_back(1.0 - d);
    }
    
    Snp& get_marker(unsigned int i) {
        return map[i];
    }
    
    string get_name(unsigned int i) const {
        return map[i].get_name();
    }
    
    double get_minor(unsigned int i) const {
        return map[i].minor();
    }

    double get_major(unsigned int i) const {
        return map[i].major();
    }
    
    double get_prob(unsigned int i, enum phased_trait pt) const {
        return map[i].get_prob(pt);
    }
    
    unsigned int get_lodscore_count() const {
        return partial_theta_count;
    }
    
    //double get_theta(unsigned int i) const ;
    //double get_inversetheta(unsigned int i) const ;
    
    inline double get_theta(unsigned int i) const {
        return thetas[i];
    }

    inline double get_inversetheta(unsigned int i) const {
        return inversethetas[i];
    }
    
    double get_theta_log(unsigned int i) const ;
    double get_inversetheta_log(unsigned int i) const ;
    
    // only work if you are using "-n 1" on the command line
    double get_theta_halfway(unsigned int i) const { return get_theta_partial(i, 1); }
    double get_inversetheta_halfway(unsigned int i) const { return get_inversetheta_partial(i, 1); }
    
    double get_genetic_position(unsigned int index, unsigned int offset) const;
    double get_theta_partial(unsigned int index, unsigned int offset) const ;
    double get_inversetheta_partial(unsigned int index, unsigned int offset) const ;
    double get_theta_partial_raw(unsigned int index) const { return partial_thetas[index]; }
    
    unsigned int num_markers() const { 
        return map.size();
    }
    
    bool sanity_check();
	string debug_string();
    
    void set_temperature(double t);
};

#endif

