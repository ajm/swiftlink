using namespace std;

#include <cmath>

#include "logarithms.h"

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
    
    return log(exp(b - a) + 1) + a;
}

double log_product(double a, double b) {
    return ((a == LOG_ZERO) or (b == LOG_ZERO)) ? LOG_ZERO : a + b;
}
