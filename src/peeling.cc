using namespace std;

#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstdio>

#include "peeling.h"

// annoyingly I cannot use random_shuffle, which might be safer (?)
PeelOperation PedigreePeeler::get_random_operation(
    vector<PeelOperation>::iterator start, vector<PeelOperation>::iterator end) {
    
    return *(start + (rand() % (end - start)));
}

// XXX add any heuristics related to peel operation selection here
// based on the last item to be added to 'peelorder' vector
// i)  if one child peeled up, then peel the rest?
// ii) if one parent peeled down, then peel the other one?
//          <-- this will probably happen automatically
PeelOperation PedigreePeeler::get_best_operation_heuristic(
    vector<PeelOperation>::iterator start, vector<PeelOperation>::iterator end) {
    return *start;
}

// XXX see Thomas'86, need to consider number of alleles
// to work out the number of computations at each step
// not just the size of the cutset
// number of alleles depends on peeling descent graph (4) or 
// genotypes during the L-sampler (3)
PeelOperation PedigreePeeler::get_best_operation(
    vector<PeelOperation>::iterator start, vector<PeelOperation>::iterator end) {

    pair<vector<PeelOperation>::iterator, PeelOperation>::iterator> bounds;
    vector<PeelOperation> v(start, end);

    sort(v.begin(), v.end());

    bounds = equal_range(v.begin(), v.end(), v[0].cutset_size());

    return get_best_operation_heuristic(bounds.first, bounds.second);
}

// create a list of all the possible peels from peripheral families
// the actual details of what a peel is, is coded in the person class
vector<PeelOperation> PedigreePeeler::all_possible_peels(unsigned int* unpeeled) {
    vector<PeelOperation> tmp;
    
    for(unsigned int i = 0; i < ped->size(); ++i) {
        PeelOperation p;

        if(not peeled[i]) {
            unpeeled++;
            
            if((*ped)[i]->peel_operation(&p)) {
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

        tmp = all_possible_peels(&unpeeled);
        
        // everyone is peeled, so we are done
        if(unpeeled == 0)
            break;

        p = get_best_operation(tmp.begin(), tmp.end());
        
        peelorder.push_back(p);
        peeled[p.get_pivot()] = true;

        tmp.clear();
    }
}

void PedigreePeeler::print() {
    
    for(unsigned int i = 0; i < peelorder.size(); ++i) {
        printf("%d: ", i);
        peelorder[i].print(); // XXX this is not ideal
    }
}

bool PedigreePeeler::peel(double *likelihood) {
    // TODO
    return false;
}

