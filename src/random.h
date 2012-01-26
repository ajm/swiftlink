#ifndef LKG_RANDOM_H_
#define LKG_RANDOM_H_


#ifdef __cplusplus
extern "C" {
#endif
    
    void seed_random(unsigned int seed);
    void destroy_random();
    double get_random();
    int get_random_int(int limit);

#ifdef __cplusplus
}
#endif

#endif
