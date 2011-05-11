#ifndef LKG_TRAIT_H_
#define LKG_TRAIT_H_

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

typedef enum unphased_trait     unphased_trait_t;
typedef enum phased_trait       phased_trait_t;

#endif

