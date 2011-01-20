#ifndef LKG_DISEASEMODEL_H_
#define LKG_DISEASEMODEL_H_

class DiseaseModel {
	double frequency;
	double penetrance[3];
	bool sexlinked;
	
 public :
	DiseaseModel() 
		: frequency(0.001), sexlinked(false) {
		penetrance[0] = 0.0;
		penetrance[1] = 0.0;
		penetrance[2] = 0.0;
	}

	double operator[](int i) {
		return penetrance[i];
	}
	
	void set_freq(double f) { frequency = f; }
	void set_penetrance(double p, const int i) { penetrance[i] = p; }
	void set_sexlinked(bool s) { sexlinked = s; }

	double get_freq() const { return frequency; }
	double get_penetrance(const int i) const { return penetrance[i]; }
	bool is_sexlinked() { return sexlinked; }

	void print() {
		printf("DiseaseModel:\n");
		printf("\tsex linked: %s\n", 
			sexlinked ? "true" : "false");
		printf("\tdisease freq: %f\n", 
			frequency);
		printf("\tpenetrance: %f, %f, %f\n", 
			penetrance[0], penetrance[1], penetrance[2]);
		printf("\n");
	}
};

#endif

