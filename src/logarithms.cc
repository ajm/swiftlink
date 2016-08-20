#include <cmath>

#include "logarithms.h"

using namespace std;


/*
double log_sum(double a, double b) {
    return log(exp(b - a) + 1) + a;
}
*/

double log_sum(double a, double b) {
    if(a == LOG_ZERO)
        return b;
    
    if(b == LOG_ZERO)
        return a;
    
    //if(b < a)
        return log(exp(b - a) + 1) + a;
    
    //return log(exp(a - b) + 1) + b;
}

double log_product(double a, double b) {
    return ((a == LOG_ZERO) or (b == LOG_ZERO)) ? LOG_ZERO : a + b;
}

double log_mean(double a, double b) {
    return log_sum(a, b) - log(2.0);
}

