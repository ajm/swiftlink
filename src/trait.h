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

enum trait {
    TRAIT_U,
    TRAIT_A
};

#define TRAIT_11 TRAIT_UU
#define TRAIT_12 TRAIT_UA
#define TRAIT_21 TRAIT_AU
#define TRAIT_22 TRAIT_AA

#define TRAIT_1 TRAIT_U
#define TRAIT_2 TRAIT_A

typedef enum unphased_trait     unphased_trait_t;
typedef enum phased_trait       phased_trait_t;

#endif

