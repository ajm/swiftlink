#ifndef LKG_DISEASEMODEL_H_
#define LKG_DISEASEMODEL_H_

#include "trait.h"
#include "person.h" // for enum affection


enum simple_disease_model {
    SIMPLE_AUTOSOMAL_RECESSIVE,
    SIMPLE_AUTOSOMAL_DOMINANT
    // XXX more...
};

class DiseaseModel {
	double frequency;
	double penetrance[3];
    double apriori_prob[3][3];
	double penetrance_prob[3][3];
	bool sexlinked;

    double get_prob(double prob[3][3], enum affection a, enum unphased_trait t);    
    void set_autosomal_recessive();
	
 public :
	DiseaseModel() 
		: frequency(0.001), sexlinked(false) {
		penetrance[0] = penetrance[1] = penetrance[2] = 0.0;
	}
    ~DiseaseModel() {}

	double operator[](int i) {
		return penetrance[i];
	}
	
	void set_freq(double f) { frequency = f; }
	void set_penetrance(double p, const enum unphased_trait t) { penetrance[t] = p; }
	void set_sexlinked(bool s) { sexlinked = s; }

	double get_freq() const { return frequency; }
	double get_penetrance(const enum unphased_trait t) const { return penetrance[t]; }
	bool is_sexlinked() { return sexlinked; }
    
    double get_penetrance_prob(enum affection a, enum unphased_trait t);
    double get_apriori_prob(enum affection a, enum unphased_trait t);
    void finish_init();

    void set(enum simple_disease_model d);
    
	void print();
};

#endif

