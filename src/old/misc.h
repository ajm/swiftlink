#ifndef LKG_MISC_H_
#define LKG_MISC_H_

#include <cstdio>
#include <cstdarg>

// used in a few places following nomenclature of
// Corman, Leiserson, Rivest & Stein 2nd Ed. 
// depth-first search
enum {
	WHITE,
	GREY,
	BLACK
};

typedef unsigned int uint;

const unsigned int DEBUG_FP_PRECISION = 3;

#define ERRORS(fmt) \
    fprintf(stderr, "error : " fmt "\n");

#define ERRORF(fmt, args...) \
    fprintf(stderr, "error : " fmt "\n", ##args);

#define DEBUGS(fmt) \
    fprintf(stderr, "error : " fmt " (%s:%d)\n", __FILE__, __LINE__);

#define DEBUGF(fmt, args...) \
    fprintf(stderr, "error : " fmt " (%s:%d)\n", ##args, __FILE__, __LINE__);

#endif

