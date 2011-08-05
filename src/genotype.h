#ifndef LKG_GENOTYPE_H_
#define LKG_GENOTYPE_H_

using namespace std;

#include <string>


enum phased_genotype {
    AA = 8,
    AB = 4,
    BA = 2,
    BB = 1,
    UN = 0
};

enum unphased_genotype {
	UNTYPED,
	HETERO,
	HOMOZ_A,
	HOMOZ_B
};

typedef enum phased_genotype    phased_genotype_t;
typedef enum unphased_genotype  unphased_genotype_t;

bool genotype_compatible(unphased_genotype_t mother, unphased_genotype_t father, 
                         unphased_genotype_t child);
bool genotype_untyped(unphased_genotype_t g);
bool genotype_homoz(unphased_genotype_t g);
bool genotype_hetero(unphased_genotype_t g);

string genotype_string(enum phased_genotype g);
string genotype_string(enum unphased_genotype g);

// add to bitmap        |=
// toggle from bitmap   ^=
// set only             &=

// TODO replace these with some inline stuff?

#define genotype_set_all(g)         ((g)  = (AA | AB | BA | BB))
#define genotype_set_homozA(g)		((g) &= AA)
#define genotype_set_homozB(g)		((g) &= BB)
#define genotype_set_hetero(g)      ((g) &= (AB | BA))
#define genotype_bad(g)			((g) == 0)
#define genotype_possible(g,v)  ((g) & (v))
#define genotype_remove(g,v)    ((g) ^= (v))
#define genotype_only(g,v)      (((g) ^ (v)) == 0)
#define genotype_single(g)		(((g) & ((g) - 1)) == 0)
#define maternal_A(g)			((g) & (AB | AA))
#define maternal_B(g)			((g) & (BA | BB))
#define paternal_A(g)			((g) & (BA | AA))
#define paternal_B(g)			((g) & (AB | BB))

#endif

