using namespace std;

#include <cstdlib>

#include "types.h"


inline void seed_random(unsigned int seed) {
    srandom(seed);
}

inline double get_random() {
    return random() / DBL_RAND_MAX;
}

inline int get_random(int limit) {
    return get_random() * limit;
}

