using namespace std;

#include <vector>
#include <algorithm>

#include "peeling.h"

// XXX  see Thomas'86, need to consider number of alleles
// to work out the number of computations at each step
// not the size of the cutset
PeelOperation PedigreePeeler::get_minimum_cutset(vector<PeelOperation> peels) {
    int tmp, tmp2;
    
    tmp = peels[0].cutset_size();
    
    for(unsigned int i = 1; i < peels.size(); ++i) {
        tmp2 = peels[i].cutset_size();
        
        if(tmp2 < tmp) {
            tmp = tmp2;
        }
    }
    
    vector<PeelOperation> nextpeels;
    
    for(unsigned int i = 0; i < peels.size(); ++i) {
        if(tmp == peels[i].cutset_size()) {
            nextpeels.push_back(peels[i]);
        }
    }
    
    // XXX reasons to prefer one operation over another?
    // i) if one child peeled up, then peel the rest?
    // ii) if one parent peeled down, then peel the other one? <-- this will probably happen automatically
    // iii) ???
    return nextpeels[0];
}

vector<PeelOperation> PedigreePeeler::get_possible_peels(unsigned int* unpeeled) {
    vector<PeelOperation> tmp;
    
    for(unsigned int i = 0; i < ped.size(); ++i) {
        PeelOperation p;

        if(not peeled[i]) {
            unpeeled++;
            
            if(ped[i]->peel_operation(&p)) {
                tmp.push_back(p);
            }
        }
    }
    
    return tmp;
}

bool PedigreePeeler::build_peel_order() {
    PeelOperation p;
    vector<PeelOperation> tmp;
    int unpeeled;
    
    fill(peeled->start(), peeled->end(), false);
    
    while(true) {
        unpeeled = 0;

        tmp = get_possible_peels(&unpeeled);
        
        // everyone is peeled, so we are done
        if(unpeeled == 0)
            break;

        p = get_minimum_cutset(tmp);
        
        peelorder.push_back(p);
        peeled[p.get_pivot()] = true;

        tmp.clear();
    }
}

