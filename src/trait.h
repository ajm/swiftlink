#ifndef LKG_TRAIT_H_
#define LKG_TRAIT_H_

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

#endif
