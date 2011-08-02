#include <cstdlib>
#include "genotype.h"

// both parents hetero - child anything
// one parent homoz - child cannot be other homoz
// both parents same homoz - child can only be homoz
// both parents different homoz - child can only be hetero
bool genotype_compatible(unphased_genotype_t mother, 
		 				 unphased_genotype_t father, 
						 unphased_genotype_t child) {
	
	switch ( child ) {
		case UNTYPED :
			return 1;
		case HOMOZ_A :
			return (mother != HOMOZ_B) && (father != HOMOZ_B);
		case HETERO :
			return !((mother == HOMOZ_A) && (father == HOMOZ_A)) && \
                   !((mother == HOMOZ_B) && (father == HOMOZ_B));
		case HOMOZ_B :
			return (mother != HOMOZ_A) && (father != HOMOZ_A);
	}
	
	abort();
}

bool genotype_untyped(unphased_genotype_t g) {
	return g == UNTYPED;
}

bool genotype_homoz(unphased_genotype_t g) {
	return g == HOMOZ_A || g == HOMOZ_B;
}

bool genotype_hetero(unphased_genotype_t g) {
	return g == HETERO;
}

string genotype_string(enum phased_genotype g) {
	switch(g) {
		case UN :
			return "UN";
		case AA :
			return "AA";
		case AB :
			return "AB";
		case BA :
			return "BA";
		case BB :
			return "BB";
	}
			
	abort();
}

string genotype_string(enum unphased_genotype g) {
    switch(g) {
        case HOMOZ_AA :
            return "AA";
        case HOMOZ_BB :
            return "BB";
        case HETERO :
            return "AB";
    }
    
    abort();
}


