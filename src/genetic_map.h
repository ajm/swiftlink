#ifndef LKG_GENETICMAP_H_
#define LKG_GENETICMAP_H_

using namespace std;

#include <cstdio>
#include <cmath>
#include <string>
#include <vector>


// XXX create a super class 'marker'
// or maybe I will just ignore more polymorphic markers
class Snp {

    string name;
    double genetic_distance;
    unsigned int physical_distance;
    double major_freq;
    double minor_freq;

 public :
    Snp(string name, double genetic, unsigned int physical) 
        : name(name), 
        genetic_distance(genetic), 
        physical_distance(physical),
        major_freq(1.0), 
        minor_freq(0.0) {}

    ~Snp() {}
    
    double major() { return major_freq; }
    double minor() { return minor_freq; }
    void set_minor_freq(double m) {
        minor_freq = m;
        major_freq = 1.0 - m;
    }
    
    string get_name() { return name; }
    double get_g_distance() { return genetic_distance; }
    unsigned int get_p_distance() { return physical_distance; }

	void print() {
		printf("\t%s\t%f\t%u\t(minor=%f, major=%f)\n",
			name.c_str(), genetic_distance, physical_distance,
			minor_freq, major_freq);
	}
};

class GeneticMap {

    vector<Snp> map;
    vector<double> thetas; // both stored in log e
    vector<double> inverse_thetas;

    double haldane(double m);
    
 public :
    GeneticMap() {}
    ~GeneticMap() {}
    
	Snp& operator[](int i) {
		return map[i];
	}

    void add(Snp& s) {
        map.push_back(s);
    }
    
    void add_theta(double d) {
        thetas.push_back(log(d));
        inverse_thetas.push_back(log1p(-d));
    }
    
    Snp& get_marker(unsigned int i) { 
        return map[i];
    }
    
    double get_minor(unsigned int i) {
        return map[i].minor();
    }

    double get_major(unsigned int i) {
        return map[i].major();
    }

    double get_theta(unsigned int i) {
        return thetas[i];
    }

    double get_inverse_theta(unsigned int i) {
        return inverse_thetas[i];
    }
    
    unsigned int num_markers() { 
        return map.size();
    }

    bool sanity_check();
	void print();
    double get_theta_halfway(unsigned int i);
};

#endif

