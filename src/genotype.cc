#include <cstdio>
#include <cstdlib>

#include "genotype.h"

// both parents hetero - child anything
// one parent homoz - child cannot be other homoz
// both parents same homoz - child can only be homoz
// both parents different homoz - child can only be hetero
bool genotype_compatible(enum unphased_genotype mother, 
		 				 enum unphased_genotype father, 
						 enum unphased_genotype child,
                         enum sex child_sex,
                         bool sex_linked) {
	
    if(not sex_linked or child_sex == FEMALE) {
	    switch (child) {
		    case UNTYPED :
			    return true;
		    case HOMOZ_A :
			    return (mother != HOMOZ_B) and (father != HOMOZ_B);
            case HOMOZ_B :
			    return (mother != HOMOZ_A) and (father != HOMOZ_A);
		    case HETERO :
			    return not ((mother == HOMOZ_A) and (father == HOMOZ_A)) and \
                       not ((mother == HOMOZ_B) and (father == HOMOZ_B));
	    }
    }
    else { // male X-linked
        switch(child) {
            case UNTYPED :
                return true;
            case HETERO :
                return false;
            case HOMOZ_A :
                return mother != HOMOZ_B;
            case HOMOZ_B :
                return mother != HOMOZ_A;
        }
    }
	
	abort();
}

bool genotype_untyped(enum unphased_genotype g) {
	return g == UNTYPED;
}

bool genotype_homoz(enum unphased_genotype g) {
	return (g == HOMOZ_A) or (g == HOMOZ_B);
}

bool genotype_hetero(enum unphased_genotype g) {
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
        case HOMOZ_A :
            return "AA";
        case HOMOZ_B :
            return "BB";
        case HETERO :
            return "AB";
        case UNTYPED :
            return "UN";
    }
    
    abort();
}

enum phased_genotype genotype_from_trait(int i) {
    
    switch(i) {
        case TRAIT_UU: return AA;
        case TRAIT_AU: return BA;
        case TRAIT_UA: return AB;
        case TRAIT_AA: return BB;
    }
    
    fprintf(stderr, "error in genotype_from_trait, %d\n", i);
    abort();
}

