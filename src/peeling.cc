using namespace std;

#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstdio>

#include "peeling.h"

// annoyingly I cannot use random_shuffle, which might be safer (?)
PeelOperation PedigreePeeler::get_random_operation(vector<PeelOperation>& v) {
    return v[rand() % v.size()];
}

// XXX add any heuristics related to peel operation selection here
// based on the last item to be added to 'peelorder' vector
// i)  if one child peeled up, then peel the rest?
// ii) if one parent peeled down, then peel the other one?
//          <-- this will probably happen automatically
PeelOperation PedigreePeeler::get_best_operation_heuristic(vector<PeelOperation>& v) {
    return v[0];
}

// XXX see Thomas'86, need to consider number of alleles
// to work out the number of computations at each step
// not just the size of the cutset
// number of alleles depends on peeling descent graph (4) or 
// genotypes during the L-sampler (3)
PeelOperation PedigreePeeler::get_best_operation(vector<PeelOperation>& v) {    

    sort(v.begin(), v.end());

    vector<PeelOperation>::iterator it = v.begin();
    unsigned int cs_size = v[0].get_cutset_size();
    while(it != v.end()) {
        if(it->get_cutset_size() != cs_size) {
            break;
        }
        it++;
    }

    vector<PeelOperation> tmp(v.begin(), it);

    return get_best_operation_heuristic(tmp);
}

// create a list of all the possible peels from peripheral families
// the actual details of what a peel is, is coded in the person class
vector<PeelOperation> PedigreePeeler::all_possible_peels(int* unpeeled) {
    vector<PeelOperation> tmp;
    Person* per;

    for(unsigned int i = 0; i < ped->num_members(); ++i) {
        PeelOperation p;

        if(not state.is_peeled(i)) {
            (*unpeeled)++;
            
            per = ped->get_by_index(i);

            if(per->peel_operation(&p, state)) {
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

    while(true) {
        unpeeled = 0;

        tmp = all_possible_peels(&unpeeled);
        
        // everyone is peeled, so we are done
        if(unpeeled == 0)
            break;

        p = get_best_operation(tmp);
/*
        // debug
        printf("step %d\n", peelorder.size());
        printf("selected: ");
        p.print();
        printf("\n");
        for(unsigned int i = 0; i < tmp.size(); ++i) {
            tmp[i].print();
        }
        printf("\n\n");
        // debug
*/
        peelorder.push_back(p);
        state.set_peeled(p.get_pivot());

        tmp.clear();
    }
}

void PedigreePeeler::print() {
    
    for(unsigned int i = 0; i < peelorder.size(); ++i) {
        printf("%d: ", i);
        peelorder[i].print(); // XXX this is not ideal
    }
}
/*
bool PedigreePeeler::peel(double *likelihood) {
    // TODO
    return false;
}
*/

