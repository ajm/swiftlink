#ifndef LKG_TYPES_H_
#define LKG_TYPES_H_

#include <cstdlib>

#include "genotype.h"
#include "trait.h"


const unsigned int UNKNOWN_PARENT = ~0u;
const unsigned int UNKNOWN_ID = UNKNOWN_PARENT;
const unsigned int DEBUG_FP_PRECISION = 3;
const double DBL_RAND_MAX = static_cast<double>(RAND_MAX);

enum parentage { 
    MATERNAL,
    PATERNAL
};

enum sex {
    UNSEXED,
    MALE,
    FEMALE
};

enum affection {
    UNKNOWN_AFFECTION,
    UNAFFECTED,
    AFFECTED
};

enum simple_disease_model {
    AUTOSOMAL_RECESSIVE,
    AUTOSOMAL_DOMINANT
};

// used in a few places following nomenclature of
// Corman, Leiserson, Rivest & Stein 2nd Ed. depth-first search
enum {
    WHITE,
    GREY,
    BLACK
};

typedef int meiosis_indicator_t;
typedef enum parentage allele_t;

string gender_str(enum sex s);
string affection_str(enum affection a);
string parent_str(enum parentage p);


#endif

