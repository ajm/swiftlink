#include <string>
#include <cstdlib>

#include "types.h"

using namespace std;


string gender_str(enum sex s) {
    
	switch(s) {
		case MALE :
			return "male";
		case FEMALE :
			return "female";
		case UNSEXED :
            return "unspecified";
	}
    
	abort();
}

string affection_str(enum affection a) {
    
	switch(a) {
		case UNAFFECTED :
			return "unaffected";
		case AFFECTED :
			return "affected";
		case UNKNOWN_AFFECTION :
            return "unspecified";
	}
    
	abort();
}

string parent_str(enum parentage p) {
    
    switch(p) {
        case MATERNAL :
            return "maternal";
        case PATERNAL :
            return "paternal";
        case NONE :
            return "none";
    }
    
    abort();
}
