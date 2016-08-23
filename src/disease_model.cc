#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
//#include <iomanip>

#include "types.h"
#include "disease_model.h"

using namespace std;


string DiseaseModel::debug_string() {
    stringstream ss;
    
    ss << "DiseaseModel:" << "\n";
    ss << "    sex linked: " << (sexlinked ? "true" : "false") << "\n";
    ss << "    disease freq: " << frequency << "\n";
    ss << "    penetrance: " << penetrance[0] << ", " 
                           << penetrance[1] << ", " 
                           << penetrance[2] << "\n";
    ss << "\n";

    return ss.str();
    
    // I hate not doing this inline, like a fmt string...
    // I know about setprecision(n) inline, but that set it permanently
    // which is 'surprising'
    ss.precision(DEBUG_FP_PRECISION);
    
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {        
            ss  << "\tapriori_prob[" << affection_str(static_cast<enum affection>(i)) << "][" \
                << trait_str(static_cast<enum unphased_trait>(j)) << "] = " \
                << scientific << apriori_prob[i][j] << "\n";
        }
    }
    
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {        
            ss  << "\tpenetrance_prob[" \
                << affection_str(static_cast<enum affection>(i)) << "][" \
                << trait_str(static_cast<enum unphased_trait>(j)) << "] = " \
                << scientific << penetrance_prob[i][j] << "\n";
        }
    }
    
    ss << "\n";
    
    return ss.str();
}

// mostly for testing
void DiseaseModel::set_autosomal_recessive() {
    
    set_freq(1 / 100000.0);
    set_sexlinked(false);
	set_penetrance(0.0, TRAIT_HOMO_U);
    set_penetrance(0.0, TRAIT_HETERO);
    set_penetrance(1.0, TRAIT_HOMO_A);

    finish_init();
}

void DiseaseModel::set_autosomal_dominant() {
    
    set_freq(1 / 100000.0);
    set_sexlinked(false);
	set_penetrance(0.0, TRAIT_HOMO_U);
    set_penetrance(1.0, TRAIT_HETERO);
    set_penetrance(1.0, TRAIT_HOMO_A);

    finish_init();
}

void DiseaseModel::set(enum simple_disease_model d) {
    
    switch(d) {
        case AUTOSOMAL_RECESSIVE :
            set_autosomal_recessive();
            break;
            
        case AUTOSOMAL_DOMINANT :
            set_autosomal_dominant();
            break;
            
        default :
            abort();
    }
    
    abort();
}

double DiseaseModel::get_prob(const double prob[3][3], enum affection a, enum unphased_trait t) const {
    return prob[a][t];
}

double DiseaseModel::get_penetrance_prob(enum affection a, enum unphased_trait t) const {
    return get_prob(penetrance_prob, a, t);
}

double DiseaseModel::get_apriori_prob(enum affection a, enum unphased_trait t) const {
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

double DiseaseModel::get_penetrance_prob2(enum affection a, enum unphased_trait t, enum sex s) const {

    if(sexlinked and s == MALE) {
        if(t == TRAIT_HETERO)
            return 0.0;
    }

    if(a == UNKNOWN_AFFECTION) {
        return (sexlinked and s == MALE) ? 0.5 : 0.25;
    }

    return ((a == AFFECTED) ? penetrance[t] : 1.0 - penetrance[t]);
}

double DiseaseModel::get_apriori_prob2(enum affection a, enum unphased_trait t, enum sex s) const {
    double tmp = 0.0;

    if(sexlinked and s == MALE) {
        if(t == TRAIT_HETERO)
            return 0.0;
    }

    switch(t) {
        case TRAIT_HOMO_A :
            tmp = (sexlinked and s == MALE) ? frequency : frequency * frequency;
            break;

        case TRAIT_HETERO :
            tmp = (sexlinked and s == MALE) ? 0.0 : frequency * (1.0 - frequency);
            break;

        case TRAIT_HOMO_U :
            tmp = (sexlinked and s == MALE) ? (1.0 - frequency) : (1.0 - frequency) * (1.0 - frequency);
            break;
    }

    switch(a) {
        case AFFECTED :
            return tmp * penetrance[t];

        case UNAFFECTED :
            return tmp * (1.0 - penetrance[t]);

        case UNKNOWN_AFFECTION :
            return (tmp * penetrance[t]) + (tmp * (1.0 - penetrance[t]));
    }

    abort();
}

