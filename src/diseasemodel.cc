using namespace std;

#include <cstdio>

#include "diseasemodel.h"
#include "trait.h"


void DiseaseModel::print() {
    printf("DiseaseModel:\n");
	printf("\tsex linked: %s\n", sexlinked ? "true" : "false");
	printf("\tdisease freq: %f\n", frequency);
	printf("\tpenetrance: %f, %f, %f\n", penetrance[0], penetrance[1], penetrance[2]);
	printf("\n");
}

enum unphased_trait DiseaseModel::phased2unphased(enum phased_trait t) {

    switch(t) {
        case TRAIT_UU :
            return TRAIT_HOMO_U;
        case TRAIT_AU :
        case TRAIT_UA :
            return TRAIT_HETERO;
        case TRAIT_AA :
            return TRAIT_HOMO_A;
        default :
            abort();
    }

    // to shutup compiler warnings
    return TRAIT_HOMO_U;
}

double DiseaseModel::get_prob(double& prob[3][3], enum affection a, enum phased_trait t) {
    return prob[a][phased2unphased(t)];
}

double DiseaseModel::get_penetrance_prob(enum affection a, enum phased_trait t) {
    return get_prob(&penetrance_prob, a, t);
}

double DiseaseModel::get_apriori_prob(enum affection a, enum phased_trait t) {
    return get_prob(&apriori_prob, a, t);
}

// this is essentially an expansion of the disease model to distributions
// of what the trait might be in an individual given their affection status
// this info will be used by the Person objects
void DiseaseModel::init_probs() {
	double affected_total;
    double unaffected_total;
    double unknown_total;
    
/*
    // don't actually normalise these because the penetrace is the probability
    // of developing the phenotype, given that you have the genetic trait
    // so they don't necessarily sum to 1
    affected_total = 1.0;
    unaffected_total = 1.0;
*/

    // XXX no, wrong. this is actually "given the affection status" what is the
    // trait allele, so this does need to be normalised
    affected_total =        penetrance[TRAIT_HOMO_A] + \
                     (2.0 * penetrance[TRAIT_HETERO]) + \
                            penetrance[TRAIT_HOMO_U];
    unaffected_total = 4.0 - affected_total;


	// penetrances + normalise inline
	penetrance_prob[AFFECTED][TRAIT_HOMO_A] = penetrance[TRAIT_HOMO_A] / affected_total;
	penetrance_prob[AFFECTED][TRAIT_HETERO] = penetrance[TRAIT_HETERO] / affected_total;
	penetrance_prob[AFFECTED][TRAIT_HOMO_U] = penetrance[TRAIT_HOMO_U] / affected_total;

	penetrance_prob[UNAFFECTED][TRAIT_HOMO_A] = (1.0 - penetrance[TRAIT_HOMO_A]) / unaffected_total;
	penetrance_prob[UNAFFECTED][TRAIT_HETERO] = (1.0 - penetrance[TRAIT_HETERO]) / unaffected_total;
	penetrance_prob[UNAFFECTED][TRAIT_HOMO_U] = (1.0 - penetrance[TRAIT_HOMO_U]) / unaffected_total;

	penetrance_prob[UNKNOWN_AFFECTION][TRAIT_HOMO_A] = 0.25;
	penetrance_prob[UNKNOWN_AFFECTION][TRAIT_HETERO] = 0.25;
	penetrance_prob[UNKNOWN_AFFECTION][TRAIT_HOMO_U] = 0.25;
    
    
	// apriori + normalise
	apriori_prob[AFFECTED][TRAIT_HOMO_A] =		  frequency  *		   frequency  * penetrance[TRAIT_HOMO_A];
	apriori_prob[AFFECTED][TRAIT_HETERO] =		  frequency  *	(1.0 - frequency) * penetrance[TRAIT_HETERO];
	apriori_prob[AFFECTED][TRAIT_HOMO_U] = (1.0 - frequency) *	(1.0 - frequency) * penetrance[TRAIT_HOMO_U];
	
	apriori_prob[UNAFFECTED][TRAIT_HOMO_A] = 	    frequency  *		   frequency  * (1.0 - penetrance[TRAIT_HOMO_A]);
	apriori_prob[UNAFFECTED][TRAIT_HETERO] = 	    frequency  *	(1.0 - frequency) * (1.0 - penetrance[TRAIT_HETERO]);
	apriori_prob[UNAFFECTED][TRAIT_HOMO_U] = (1.0 - frequency) *	(1.0 - frequency) * (1.0 - penetrance[TRAIT_HOMO_U]);
	
	apriori_prob[UNKNOWN_AFFECTION][TRAIT_HOMO_A] = apriori_prob[AFFECTED][TRAIT_HOMO_A] + apriori_prob[UNAFFECTED][TRAIT_HOMO_A];
	apriori_prob[UNKNOWN_AFFECTION][TRAIT_HETERO] = apriori_prob[AFFECTED][TRAIT_HETERO] + apriori_prob[UNAFFECTED][TRAIT_HETERO];
	apriori_prob[UNKNOWN_AFFECTION][TRAIT_HOMO_U] = apriori_prob[AFFECTED][TRAIT_HOMO_U] + apriori_prob[UNAFFECTED][TRAIT_HOMO_U];

    // sum
    affected_total =       apriori_prob[AFFECTED][TRAIT_HOMO_A] + \
					(2.0 * apriori_prob[AFFECTED][TRAIT_HETERO]) + \
					       apriori_prob[AFFECTED][TRAIT_HOMO_U];
	
	unaffected_total =     apriori_prob[UNAFFECTED][TRAIT_HOMO_A] + \
					(2.0 * apriori_prob[UNAFFECTED][TRAIT_HETERO]) + \
						   apriori_prob[UNAFFECTED][TRAIT_HOMO_U];

    unknown_total =        apriori_prob[UNKNOWN_AFFECTION][TRAIT_HOMO_A] + \
                    (2.0 * apriori_prob[UNKNOWN_AFFECTION][TRAIT_HETERO]) + \
                           apriori_prob[UNKNOWN_AFFECTION][TRAIT_HOMO_U];
	
    // normalise
	apriori_prob[AFFECTED][TRAIT_HOMO_A] /= affected_total;
	apriori_prob[AFFECTED][TRAIT_HETERO] /= affected_total;
	apriori_prob[AFFECTED][TRAIT_HOMO_U] /= affected_total;
	
	apriori_prob[UNAFFECTED][TRAIT_HOMO_A] /= unaffected_total;
	apriori_prob[UNAFFECTED][TRAIT_HETERO] /= unaffected_total;
	apriori_prob[UNAFFECTED][TRAIT_HOMO_U] /= unaffected_total;

    apriori_prob[UNKNOWN_AFFECTION][TRAIT_HOMO_A] /= unknown_total;
	apriori_prob[UNKNOWN_AFFECTION][TRAIT_HETERO] /= unknown_total;
	apriori_prob[UNKNOWN_AFFECTION][TRAIT_HOMO_U] /= unknown_total;
}

