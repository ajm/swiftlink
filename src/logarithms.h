#ifndef LKG_LOGARITHMS_H_
#define LKG_LOGARITHMS_H_

using namespace std;

#include <limits>

const double LOG_ILLEGAL = -numeric_limits<double>::max();
const double LOG_ZERO    = LOG_ILLEGAL;

double log_sum(double a, double b);
double log_product(double a, double b);
double log_mean(double a, double b);

#endif

