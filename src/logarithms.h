#ifndef LKG_LOGARITHMS_H_
#define LKG_LOGARITHMS_H_

#include <cmath>


double log_max(double *vals, int count) {
    double max = vals[0];
    
    for(int i = 1; i < count; ++i) {
        if(vals[i] > max) {
            max = vals[i];
        }
    }
    
    return max;
}

double log_sum(double *vals, int count) {
    double max;
    double total;
    
    max = log_max(vals, count);
    total = 0.0;
    
    for(int i = 0; i < count; ++i) {
        total += exp(vals[i] - max);
    }
    
    return max + log(total);
}

#endif

