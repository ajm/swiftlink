using namespace std;

#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
//#include <iomanip>

#include "disease_model.h"
#include "trait.h"


string DiseaseModel::debug_string() {
    stringstream ss;
    
    ss << "DiseaseModel:" << endl;
    ss << "\tsex linked: " << (sexlinked ? "true" : "false") << endl;
    ss << "\tdisease freq: " << frequency << endl;
    ss << "\tpenetrance: " << penetrance[0] << ", " 
                           << penetrance[1] << ", " 
                           << penetrance[2] << endl;
    ss << endl;
    
    // I hate not doing this inline, like a fmt string...
    // I know about setprecision(n) inline, but that set it permanently
    // which is 'surprising'
    ss.precision(3);
    
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {        
            ss << "\tapriori_prob[" << i << "][" << j << "] = " \
               << scientific << apriori_prob[i][j] << endl;
        }
    }
    
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {        
            ss << "\tpenetrance_prob[" << i << "][" << j << "] = " \
               << scientific << penetrance_prob[i][j] << endl;
        }
    }
    
    ss << endl;
    ss << "(1st index: UNKNOWN_AFFECTION, UNAFFECTED, AFFECTED)" << endl; 
    ss << "(2nd index: TRAIT_HOMO_U, TRAIT_HETERO, TRAIT_HOMO_A)" << endl;
    ss << endl;
    
    return ss.str();
}

// mostly for testing
void DiseaseModel::set_autosomal_recessive() {
    
    set_freq(1 / 1000.0); // XXX sensible value?
    set_sexlinked(false);
	set_penetrance(0.0, TRAIT_HOMO_U);
    set_penetrance(0.0, TRAIT_HETERO);
    set_penetrance(1.0, TRAIT_HOMO_A);

    finish_init();
}

void DiseaseModel::set_autosomal_dominant() {
    
    set_freq(1 / 1000.0); // XXX sensible value?
    set_sexlinked(false);
	set_penetrance(0.0, TRAIT_HOMO_U);
    set_penetrance(1.0, TRAIT_HETERO);
    set_penetrance(1.0, TRAIT_HOMO_A);

    finish_init();
}

void DiseaseModel::set(enum simple_disease_model d) {
    switch(d) {
        case SIMPLE_AUTOSOMAL_RECESSIVE :
            set_autosomal_recessive();
            break;
        case SIMPLE_AUTOSOMAL_DOMINANT :
            set_autosomal_dominant();
            break;
        default :
            abort();
    }
}

double DiseaseModel::get_prob(double prob[3][3], enum affection a, enum unphased_trait t) {
    return prob[a][t];
}

double DiseaseModel::get_penetrance_prob(enum affection a, enum unphased_trait t) {
    return get_prob(penetrance_prob, a, t);
}

double DiseaseModel::get_apriori_prob(enum affection a, enum unphased_trait t) {
    return get_prob(apriori_prob, a, t);
}

// this is essentially an expansion of the disease model to distributions
// of what the trait might be in an individual given their affection status
// this info will be used by the Person objects
void DiseaseModel::finish_init() {
	// penetrances
	penetrance_prob[AFFECTED][TRAIT_HOMO_A] = penetrance[TRAIT_HOMO_A];
	penetrance_prob[AFFECTED][TRAIT_HETERO] = penetrance[TRAIT_HETERO];
	penetrance_prob[AFFECTED][TRAIT_HOMO_U] = penetrance[TRAIT_HOMO_U];

	penetrance_prob[UNAFFECTED][TRAIT_HOMO_A] = (1.0 - penetrance[TRAIT_HOMO_A]);
	penetrance_prob[UNAFFECTED][TRAIT_HETERO] = (1.0 - penetrance[TRAIT_HETERO]);
	penetrance_prob[UNAFFECTED][TRAIT_HOMO_U] = (1.0 - penetrance[TRAIT_HOMO_U]);

	penetrance_prob[UNKNOWN_AFFECTION][TRAIT_HOMO_A] = 0.25;
	penetrance_prob[UNKNOWN_AFFECTION][TRAIT_HETERO] = 0.25;
	penetrance_prob[UNKNOWN_AFFECTION][TRAIT_HOMO_U] = 0.25;
    
    
	// apriori
	apriori_prob[AFFECTED][TRAIT_HOMO_A] =		  frequency  *		   frequency  * penetrance[TRAIT_HOMO_A];
	apriori_prob[AFFECTED][TRAIT_HETERO] =		  frequency  *	(1.0 - frequency) * penetrance[TRAIT_HETERO];
	apriori_prob[AFFECTED][TRAIT_HOMO_U] = (1.0 - frequency) *	(1.0 - frequency) * penetrance[TRAIT_HOMO_U];
	
	apriori_prob[UNAFFECTED][TRAIT_HOMO_A] = 	    frequency  *		   frequency  * (1.0 - penetrance[TRAIT_HOMO_A]);
	apriori_prob[UNAFFECTED][TRAIT_HETERO] = 	    frequency  *	(1.0 - frequency) * (1.0 - penetrance[TRAIT_HETERO]);
	apriori_prob[UNAFFECTED][TRAIT_HOMO_U] = (1.0 - frequency) *	(1.0 - frequency) * (1.0 - penetrance[TRAIT_HOMO_U]);
	
	apriori_prob[UNKNOWN_AFFECTION][TRAIT_HOMO_A] = apriori_prob[AFFECTED][TRAIT_HOMO_A] + apriori_prob[UNAFFECTED][TRAIT_HOMO_A];
	apriori_prob[UNKNOWN_AFFECTION][TRAIT_HETERO] = apriori_prob[AFFECTED][TRAIT_HETERO] + apriori_prob[UNAFFECTED][TRAIT_HETERO];
	apriori_prob[UNKNOWN_AFFECTION][TRAIT_HOMO_U] = apriori_prob[AFFECTED][TRAIT_HOMO_U] + apriori_prob[UNAFFECTED][TRAIT_HOMO_U];

}

