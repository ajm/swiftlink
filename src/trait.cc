#include <cstdlib>
#include <string>

#include "trait.h"

using namespace std;


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

string fake_trait_str(int i) {

    if(i == -1)
        return "XX";

    enum phased_trait pt = static_cast<enum phased_trait>(i);    
    
    switch(pt) {
        case TRAIT_UU:
            return "AA";
        case TRAIT_AU :
            return "BA";
        case TRAIT_UA :
            return "AB";
        case TRAIT_AA :
            return "BB";
    }
    
    abort();
}
