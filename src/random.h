#ifndef LKG_RANDOM_H_
#define LKG_RANDOM_H_

using namespace std;

void init_random();
void destroy_random();

void seed_random_explicit(string filename);
void seed_random_implicit();

double get_random();
int get_random_int(int limit);

#endif

