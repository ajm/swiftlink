#ifndef LKG_TRAIT_H_
#define LKG_TRAIT_H_

//using namespace std;

//#include <cstdlib>


enum linkagefile_loci {
    QUANTITATIVE_VARIABLE,
    AFFECTION_STATUS,
    BINARY_FACTOR,
    NUMBERED_ALLELES
};

// these need to be ordered the same as penetrance vector
enum unphased_trait {
    TRAIT_HOMO_U,   // phenocopies
    TRAIT_HETERO,
    TRAIT_HOMO_A
};

enum phased_trait {
    TRAIT_UU,
    TRAIT_AU,
    TRAIT_UA,
    TRAIT_AA
};

typedef enum linkagefile_loci   linkagefile_loci_t;
typedef enum unphased_trait     unphased_trait_t;
typedef enum phased_trait       phased_trait_t;

/*
enum unphased_trait phased2unphased(enum phased_trait t) {

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

const char *loci2str(const linkagefile_loci_t l) {
    switch (l) {
        case QUANTITATIVE_VARIABLE :    return "quantitative variable";
        case AFFECTION_STATUS :         return "affection status";
        case BINARY_FACTOR :            return "binary factor";
        case NUMBERED_ALLELES :         return "numbered allele";
        default:						break;
    }
	
    return "unknown type";
}
*/
#endif

