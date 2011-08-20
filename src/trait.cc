using namespace std;

#include <cstdlib>
#include <string>

#include "trait.h"


string trait_str(enum phased_trait pt) {
    switch(pt) {
        case TRAIT_UU :
            return "UU";
        case TRAIT_AU :
            return "AU";
        case TRAIT_UA :
            return "UA";
        case TRAIT_AA :
            return "AA";
    }
    
    abort();
}

string trait_str(enum unphased_trait pt) {
    switch(pt) {
        case TRAIT_HOMO_U :
            return "UU";
        case TRAIT_HETERO :
            return "AU";
        case TRAIT_HOMO_A :
            return "AA";
    }
    
    abort();
}
