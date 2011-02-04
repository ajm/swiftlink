#ifndef LKG_RFUNCTION_H_
#define LKG_RFUNCTION_H_

using namespace std;

#include <cmath>
#include <queue>

class Rfunction {
    
    PeelOperation& peel;
    unsigned int num_alleles;
    double* rfunc;
    
    void normalise() {} // specific dimension(s)?
    
    
 public :
    Rfunction(PeelOperation& p, unsigned int alleles) 
        : peel(p), num_alleles(alleles) {
        rfunc = new double[static_cast<int>(pow(num_alleles, peel.get_cutset_size()))];
    }
    
    ~Rfunction() {
        delete[] rfunc;
    }
    
    void evaluate() {
        queue<unsigned int> q;
        
        for(unsigned int i = 0; i < peel.get_cutset_size(); ++i) {
            q.push_back(0);
        }
        
        while(not q.empty()) {
            tmp = q.front() + 1;
            q.pop();
            
            if(tmp < num_alleles) {
                q.push_back(tmp);
            }
            
            if (q.size() == peel.get_cutset_size()) {
                // call a function to set a value
            }
            else {
                q.push_back(0);
            }
        }
    }
};

#endif

