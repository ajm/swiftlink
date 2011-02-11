#ifndef LKG_DISEASEMODEL_H_
#define LKG_DISEASEMODEL_H_

#include "trait.h"
#include "person.h" // for enum affection


class DiseaseModel {
	double frequency;
	double penetrance[3];
    double apriori_prob[3][3];
	double penetrance_prob[3][3];
	bool sexlinked;

    enum unphased_trait phased2unphased(enum phased_trait t);
    double get_prob(double& prob[3][3], enum affection a, enum phased_trait t);
	
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
	void set_penetrance(double p, const int i) { penetrance[i] = p; }
	void set_sexlinked(bool s) { sexlinked = s; }

	double get_freq() const { return frequency; }
	double get_penetrance(const int i) const { return penetrance[i]; }
	bool is_sexlinked() { return sexlinked; }
    
    double get_penetrance_prob(enum affection a, enum phased_trait t);
    double get_apriori_prob(enum affection a, enum phased_trait t);
    void init_probs();
    
	void print();
};

#endif

