using namespace std;

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include "peeling.h"
#include "peel_sequence_generator.h"


PeelOperation PeelSequenceGenerator::get_random_operation(vector<PeelOperation>& v) {
    return v[random() % v.size()];
}

// XXX add any heuristics related to peel operation selection here
// based on the last item to be added to 'peelorder' vector
// i)  if one child peeled up, then peel the rest?
// ii) if one parent peeled down, then peel the other one?
//          <-- this will probably happen automatically
PeelOperation PeelSequenceGenerator::get_best_operation_heuristic(vector<PeelOperation>& v) {

    for(unsigned int i = 0; i < v.size(); ++i) {
        if(v[i].get_type() == CHILD_PEEL) {
            return v[i];
        }
    }
    
    return v[0];
}

// XXX see Thomas'86, need to consider number of alleles
// to work out the number of computations at each step
// not just the size of the cutset
// number of alleles depends on peeling descent graph (4) or 
// genotypes during the L-sampler (3)
PeelOperation PeelSequenceGenerator::get_best_operation(vector<PeelOperation>& v) {    

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
void PeelSequenceGenerator::all_possible_peels(int& unpeeled) {
    Person* per;

    tmp.clear();

    for(unsigned int i = 0; i < ped.num_members(); ++i) {
        PeelOperation p;

        if(not state.is_peeled(i)) {
            unpeeled++;
            
            per = ped.get_by_index(i);

            if(per->peel_operation(p, state)) {
                tmp.push_back(p);
            }
        }
    }
}

bool PeelSequenceGenerator::build_peel_order() {

    while(true) {
        int unpeeled = 0;

        all_possible_peels(unpeeled); // populates 'tmp'
        
        // everyone is peeled, so we are done
        if(unpeeled == 0) {
            break;
        }

        if((unpeeled != 0) and (tmp.size() == 0)) {
            fprintf(stderr, "Error, %s failed to produce a valid peeling sequence\n", __func__ );
            return false;
        }
        
        PeelOperation p = get_best_operation(tmp);
/*
        if(p.get_type() == LAST_PEEL) {
            state.set_peeled(p.get_pivot());
            break;
        }
*/        
        peelorder.push_back(p);
        state.set_peeled(p.get_pivot());
    }

    return true;
}

vector<PeelOperation>& PeelSequenceGenerator::get_peel_order() {
    return peelorder;
}

void PeelSequenceGenerator::print() {
    
    for(unsigned int i = 0; i < peelorder.size(); ++i) {
        printf("rfunction %d: ", i);
        peelorder[i].print(); // XXX this is not ideal
    }
}

