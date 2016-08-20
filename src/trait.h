#ifndef LKG_TRAIT_H_
#define LKG_TRAIT_H_

using namespace std;

#include <string>

// these need to be ordered the same as penetrance vector
// which is the same as in the 'linkage' config file
enum unphased_trait {
    TRAIT_HOMO_U,
    TRAIT_HETERO,
    TRAIT_HOMO_A
};

enum trait {
    TRAIT_U,
    TRAIT_A
};

enum phased_trait {
    TRAIT_UU,
    TRAIT_AA,
    TRAIT_AU,
    TRAIT_UA //,
    //TRAIT_AA
};

string trait_str(enum phased_trait t);
string trait_str(enum unphased_trait t);
string fake_trait_str(int i);

#endif
