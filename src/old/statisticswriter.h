#ifndef LKG_STATSWRITER_H_
#define LKG_STATSWRITER_H_

using namespace std;

#include <cmath>
#include <fstream>
#include <algorithm>

#include "descentgraph.h"
#include "logarithms.h"

class StatisticsWriter {
 
    fstream out;
 
 public:
    StatisticsWriter() {
        out.open("statistics.txt", ios::out);
    }
    
    ~StatisticsWriter() {
        out.close();
    }
    
    void print(double* values, int count) {
        double total, mean;
        
        sort(values, values + count);
        
        total = log_sum(values, count);
        mean = total - log(count);
                
        out << mean / log(10) << '\t';
        out << values[count/2] / log(10) << '\t';
        for(int i = 0; i < count; ++i) {
            if(values[i] != LOG_ILLEGAL) {
                out << values[i] / log(10) << '\t';
                break;
            }
        }
        out << values[count-1] / log(10) << endl;
    }
};

#endif

