#ifndef LKG_DISEASEMODEL_H_
#define LKG_DISEASEMODEL_H_

using namespace std;

#include <string>

#include "types.h"


class DiseaseModel {
	double frequency;
	double penetrance[3];
    double apriori_prob[3][3];
	double penetrance_prob[3][3];
	bool sexlinked;
    
    double get_prob(const double prob[3][3], enum affection a, enum unphased_trait t) const;    
    void set_autosomal_recessive();
    void set_autosomal_dominant();
	
 public :
	DiseaseModel() : 
        frequency(0.001), 
        sexlinked(false) {
		
        penetrance[0] = penetrance[1] = penetrance[2] = 0.0;
	}

    DiseaseModel(double freq, double pen[3], bool sexlink) :
        frequency(freq),
        sexlinked(sexlink) {
  
        for(unsigned int i = 0; i < 3; ++i)
            penetrance[i] = pen[i];
    
        finish_init();    
    }
    
    ~DiseaseModel() {}
    
	void set_freq(double f) { frequency = f; }
	void set_penetrance(double p, const enum unphased_trait t) { penetrance[t] = p; }
	void set_sexlinked(bool s) { sexlinked = s; }

	double get_freq() const { return frequency; }
	double get_penetrance(const enum unphased_trait t) const { return penetrance[t]; }
	bool is_sexlinked() { return sexlinked; }
    
    double get_penetrance_prob(enum affection a, enum unphased_trait t) const;
    double get_apriori_prob(enum affection a, enum unphased_trait t) const;
    
    void set(enum simple_disease_model d);
    void finish_init();
    
	string debug_string();
};

#endif

