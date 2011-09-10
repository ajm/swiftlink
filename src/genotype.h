#ifndef LKG_GENOTYPE_H_
#define LKG_GENOTYPE_H_

using namespace std;

#include <string>

// this is kind of old-school, my initial attempt at this
// program was in C, so there are the occasional C-idioms
// in certain classes, eg: Elimination & FounderAlleleGraph

// in newer bits of the program when I need to have a phased
// genotype I just use the phased_trait enum
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

bool genotype_compatible(enum unphased_genotype mother, 
                         enum unphased_genotype father, 
                         enum unphased_genotype child);

bool genotype_untyped(enum unphased_genotype g);
bool genotype_homoz  (enum unphased_genotype g);
bool genotype_hetero (enum unphased_genotype g);

string genotype_string(enum phased_genotype g);
string genotype_string(enum unphased_genotype g);

// add to bitmap        |=
// toggle from bitmap   ^=
// set only             &=

// replace these with some inline stuff?

#define genotype_set_all(g)     ((g)  = (AA | AB | BA | BB))
#define genotype_set_homozA(g)	((g) &= AA)
#define genotype_set_homozB(g)	((g) &= BB)
#define genotype_set_hetero(g)  ((g) &= (AB | BA))
#define genotype_bad(g)         ((g) == 0)
#define genotype_possible(g,v)  ((g) & (v))
#define genotype_remove(g,v)    ((g) ^= (v))
#define genotype_only(g,v)      (((g) ^ (v)) == 0)
#define genotype_single(g)		(((g) & ((g) - 1)) == 0)
#define maternal_A(g)			((g) & (AB | AA))
#define maternal_B(g)			((g) & (BA | BB))
#define paternal_A(g)			((g) & (BA | AA))
#define paternal_B(g)			((g) & (AB | BB))

#endif

