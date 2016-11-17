#ifndef LKG_OMPFACADE_H_
#define LKG_OMPFACADE_H_

#include <cstdio>
#include <cstdlib>

#if defined(_OPENMP)
    #include <omp.h>
#endif

inline int get_max_threads() {
    #if defined(_OPENMP)
    return omp_get_max_threads();
    #else
    return 1;
    #endif
}

inline void set_num_threads(int t) {
    #if defined(_OPENMP)
    omp_set_num_threads(t);
    #else
    if(t != 1) {
        fprintf(stderr, "WARNING: this binary was compiled without OpenMP support, so does not support multithreading\n");
    }
    #endif
}

inline int get_thread_num() {
    #if defined(_OPENMP)
    return omp_get_thread_num();
    #else
    return 0;
    #endif
}

inline double get_wtime() {
    #if defined(_OPENMP)
    return omp_get_wtime();
    #else
    fprintf(stderr, "get_wtime() should not be called without multiprocessing support\n");
    abort();
    #endif
}

#endif

