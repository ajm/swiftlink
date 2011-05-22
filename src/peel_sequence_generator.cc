using namespace std;

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include "peeling.h"
#include "peel_sequence_generator.h"
#include "person.h"


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
    
    
    // XXX
    printf("\n");
    for(unsigned i = 0; i < v.size(); ++i) {
        printf("\t");
        v[i].print();
    }
    

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

bool PeelSequenceGenerator::creates_simple_peel_sequence(PeelOperation& po) {
    
    //printf("*testing: ");
    //po.print();
    
    // this is the first thing to be peeled
    if(peelorder.size() == 0)
        return true;
    
    PeelOperation& prev = peelorder.back();
    
    // in the case of CHILD_PEEL, the cutsets are the same
    if(prev.get_cutset() == po.get_cutset()) {
        return true;
    }
    
    // the cutset of prev has a member removed
    // CHILD_PEEL, PARENT_PEEL, LAST_PEEL this is the 'pivot'
    // PARENT_PEEL this is the two parents of the 'pivot'
    
    
    // CUTSETS ARE IDENTICAL _OR_
    // EVERY NODE IN PEELSET IS EITHER IN THE PREVIOUS
    // CUTSET OR ELSE IS A FOUNDER OR A LEAF
    // (THEN THE FIRST TEST CAN BE REMOVED ON LINE 70)
    
    vector<unsigned>& new_peelset = po.get_peelset();
    vector<unsigned>& old_cutset = prev.get_cutset();
    
    for(unsigned i = 0; i < new_peelset.size(); ++i) {
        if(find(old_cutset.begin(), old_cutset.end(), new_peelset[i]) == old_cutset.end()) {
            Person* per = ped.get_by_index(new_peelset[i]);
            if(not (per->isfounder() or per->isleaf())) {
                return false;
            }
        }
    }
        
    return true;
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
                if(creates_simple_peel_sequence(p)) {
                    tmp.push_back(p);
                }
            }
        }
    }
}

void PeelSequenceGenerator::build_peel_order() {

    while(true) {
        int unpeeled = 0;

        all_possible_peels(unpeeled); // populates 'tmp'
        
        // everyone is peeled, so we are done
        if(unpeeled == 0) {
            break;
        }

        if((unpeeled != 0) and (tmp.size() == 0)) {
            fprintf(stderr, "Error, %s failed to produce a valid peeling sequence\n", __func__ );
            abort();
        }
        
        PeelOperation p = get_best_operation(tmp);
/*
        if(p.get_type() == LAST_PEEL) {
            state.set_peeled(p.get_pivot());
            break;
        }
*/
        // XXX
        printf("\nselected: ");
        p.print();
        printf("\n\n");

        peelorder.push_back(p);
        
        state.toggle_peel_operation(p);
        state.print();
    }
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

