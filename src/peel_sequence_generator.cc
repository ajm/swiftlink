using namespace std;

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

#include "peeling.h"
#include "peel_sequence_generator.h"
#include "person.h"
#include "random.h"


PeelOperation PeelSequenceGenerator::get_random_operation(vector<PeelOperation>& v) {
    return v[get_random_int(v.size())];
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
    
    for(unsigned i = 0; i < ped->num_members(); ++i) {
        //PeelOperation p;
        
        if(not state.is_peeled(i)) {
            unpeeled++;
            
            per = ped->get_by_index(i);
            
            PeelOperation p = per->peel_operation(state);

            if(p.get_type() != NULL_PEEL)
                tmp.push_back(p);

            //if(per->peel_operation(p, state)) {
            //    tmp.push_back(p);
            //}
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

        find_previous_functions(p);
        
        //bruteforce_assignments(p);

        peelorder.push_back(p);
        
        state.toggle_peel_operation(p);
        
        printf("%s\n", p.debug_string().c_str());
    }
    
    /*
    // use greedy solution as input to random down hill search
    vector<int> current;
    for(unsigned i = 0; i < peelorder.size(); ++i) {
        current.push_back(peelorder[i].get_peelnode());
    }
    
    printf("peel cost start: %d\n", calculate_cost(current));
    
    int swap0, swap1, tmp, new_cost;
    int cost = calculate_cost(current);
    
    swap0 = swap1 = 0;
    
    for(unsigned int i = 0; i < 100000; ++i) {
        
        do {
            swap0 = get_random_int(current.size());
            swap1 = get_random_int(current.size());
        
        } while(swap0 == swap1);
        
        // swap random peels
        tmp = current[swap0];
        current[swap0] = current[swap1];
        current[swap1] = tmp;
        
        new_cost = calculate_cost(current);
        
        //printf("iteration %d: %d (%d, %d)\n", i, new_cost, swap0, swap1);
        
        if((new_cost != -1) and (new_cost < cost)) {
            printf("iteration %d: %d\n", i, new_cost);
        }
        
        // if better, store result
        if((new_cost != -1) and (new_cost <= cost)) {
            cost = new_cost;
            continue;
        }
        
        // not better swap back
        tmp = current[swap0];
        current[swap0] = current[swap1];
        current[swap1] = tmp;
    }
    
    printf("peel cost end: %d\n", calculate_cost(current));
    
    
    rebuild_peel_order(current);
    
    //exit(-1);
    */
    
    /*
    vector<int> v;
    
    v.push_back(10);
    v.push_back(11);
    v.push_back(1);
    v.push_back(13);
    v.push_back(0);
    v.push_back(16);
    v.push_back(12);
    v.push_back(4);
    v.push_back(14);
    v.push_back(6);
    v.push_back(27);
    v.push_back(17);
    v.push_back(18);
    v.push_back(25);
    v.push_back(2);
    v.push_back(26);
    v.push_back(19);
    v.push_back(7);
    v.push_back(22);
    v.push_back(24);
    v.push_back(3);
    v.push_back(23);
    v.push_back(20);
    v.push_back(5);
    v.push_back(21);
    v.push_back(15);
    v.push_back(9);
    v.push_back(8);
    
    rebuild_peel_order(v);
    */
    /*
    vector<int> v;
    
    v.push_back(12);
    v.push_back(6);
    v.push_back(0);
    v.push_back(5);
    v.push_back(7);
    v.push_back(9);
    v.push_back(1);
    v.push_back(8);
    v.push_back(11);
    v.push_back(4);
    v.push_back(10);
    v.push_back(2);
    v.push_back(3);
    
    rebuild_peel_order(v);
    */
    
    for(unsigned int i = 0; i < peelorder.size(); ++i) {
        bruteforce_assignments(peelorder[i]);
    }
}

vector<PeelOperation>& PeelSequenceGenerator::get_peel_order() {
    return peelorder;
}

unsigned PeelSequenceGenerator::score_peel_sequence() {
    unsigned total = 0;
    
    for(unsigned i = 0; i < peelorder.size(); ++i) {
        total += pow(4.0, static_cast<double>(peelorder[i].get_cutset_size()));
    }
    
    return total;
}

string PeelSequenceGenerator::debug_string() {
    stringstream ss;
    
    for(unsigned i = 0; i < peelorder.size(); ++i) {
        ss << peelorder[i].debug_string() << "\n";
    }
    
    return ss.str();
}

void PeelSequenceGenerator::find_previous_functions(PeelOperation& op) {
    if(op.get_type() == CHILD_PEEL) {
        find_child_functions(op);
    }
    else {
        find_generic_functions(op);
    }
}

void PeelSequenceGenerator::find_generic_functions(PeelOperation& op) {
    vector<unsigned> tmp;
    tmp.push_back(op.get_peelnode());
    
    op.set_previous_operation(find_function_containing(tmp));
    op.set_previous_operation(find_function_containing(tmp));
    
    // how can i tell which is anc and which is des
}

void PeelSequenceGenerator::find_child_functions(PeelOperation& op) {
    Person* p = ped->get_by_index(op.get_peelnode());
    vector<unsigned> tmp;
    tmp.push_back(p->get_maternalid());
    tmp.push_back(p->get_paternalid());
    
    // ancestors
    op.set_previous_operation(find_function_containing(tmp));
    
    if(p->isleaf())
        return;
    
    tmp.clear();
    tmp.push_back(p->get_internalid());
    
    // descendents
    op.set_previous_operation(find_function_containing(tmp));
}

int PeelSequenceGenerator::find_function_containing(vector<unsigned>& nodes) {
    for(int i = 0; i < int(peelorder.size()); ++i) {
        if(peelorder[i].is_used())
            continue;
        
        if(peelorder[i].contains_cutnodes(nodes)) {
            peelorder[i].set_used();
            return i;
        }
    }
    
    return -1;
}

void PeelSequenceGenerator::bruteforce_assignments(PeelOperation& op) {
    int ndim = op.get_cutset_size();
    int total = pow(4.0, ndim);
    int offset, index;
    
    vector<unsigned int>& cutset = op.get_cutset();
    
    vector<vector<int> > assigns(total, vector<int>(ped->num_members(), -1));
    
    for(int ind = 0; ind < total; ++ind) {
        index = ind;
        
        for(int i = ndim - 1; i > -1; --i) {
            offset = 1 << (i * 2);
            assigns[ind][cutset[i]] = index / offset;
            index %= offset;
        }
        
        /*
        for(int i = 0; i < ndim; ++i) {
            fprintf(stderr, "%d ", assigns[ind][cutset[i]]);
        }
        fprintf(stderr, "\n");
        */
    }
    
    //fprintf(stderr, "===\n");
    
    op.set_index_values(assigns);
}

int PeelSequenceGenerator::calculate_cost(vector<int>& seq) {
    PeelingState ps(ped);
    int cost = 0;
    
    for(unsigned i = 0; i < seq.size(); ++i) {
        Person* per = ped->get_by_index(seq[i]);
        PeelOperation p = per->peel_operation(ps);
        
        if(p.get_type() == NULL_PEEL) {
            return -1;
        }
        
        ps.set_peeled(seq[i]);
        
        cost += pow(4.0, int(p.get_cutset_size()));
        
        //printf("%d %d\n", seq[i], p->get_cutset_size(ps));
    }
    
    return cost;
}

void PeelSequenceGenerator::rebuild_peel_order(vector<int>& seq) {
    
    peelorder.clear();
    state.reset();

    for(unsigned i = 0; i < seq.size(); ++i) {
        Person* per = ped->get_by_index(seq[i]);
        PeelOperation p = per->peel_operation(state);
        
        if(p.get_type() == NULL_PEEL) {
            printf("bad rebuild: %s\n", p.debug_string().c_str());
            abort();
        }
        
        find_previous_functions(p);

        printf("rebuild: %s\n", p.debug_string().c_str());
        
        state.toggle_peel_operation(p);
        
        peelorder.push_back(p);
    }
}

