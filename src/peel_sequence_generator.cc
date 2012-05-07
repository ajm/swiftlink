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
#include "progress.h"
#include "simple_parser.h"
#include "elimination.h"
#include "genetic_map.h"


PeelOperation PeelSequenceGenerator::get_random_operation(vector<PeelOperation>& v) {
    return v[get_random_int(v.size())];
}

PeelOperation PeelSequenceGenerator::get_best_operation_heuristic(vector<PeelOperation>& v) {
    
    for(unsigned int i = 0; i < v.size(); ++i) {
        if(v[i].get_type() == CHILD_PEEL) {
            return v[i];
        }
    }
    
    return v[0];
}

PeelOperation PeelSequenceGenerator::get_best_operation(vector<PeelOperation>& v) {    

    sort(v.begin(), v.end());

    vector<PeelOperation>::iterator it = v.begin();
    unsigned int cs_size = v[0].get_cutset_size();
    while(it != v.end()) {
        if(it->get_cutset_size() != cs_size) {
            break;
        }
        ++it;
    }

    vector<PeelOperation> tmp(v.begin(), it);

    //return get_best_operation_heuristic(tmp);
    return get_random_operation(tmp);
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

    int start_index = 0;
    bool use_random = false;

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
        
        //PeelOperation p = get_best_operation(tmp);
        
        
        PeelOperation p = use_random ? get_random_operation(tmp) : get_best_operation(tmp);
        
        if((not use_random) and (p.get_cutset_size() > 2)) {
            start_index = peelorder.size();
            use_random = true;
            printf("[switching to random]\n");
            continue;
        }
        
        
        find_previous_functions(p);
        
        //bruteforce_assignments(p);
        
        peelorder.push_back(p);
        
        state.toggle_peel_operation(p);
        
        printf("%s\n", p.debug_string().c_str());
    }
        
    
    int iterations = 1000000;
    Progress p("Peel Sequence:", iterations);
    
    // section 6.3, peeling sequence generation
    // use greedy solution as input to random down hill search
    vector<unsigned int> current;
    for(unsigned i = 0; i < peelorder.size(); ++i) {
        current.push_back(peelorder[i].get_peelnode());
    }
    
    
    //printf("peel cost start: %d\n", calculate_cost(current));
    
    
    int swap0, swap1, tmp, new_cost;
    int cost = calculate_cost(current);
    int iter = -1;
    //int best_cost = cost;
    
    swap0 = swap1 = 0;
    
    int total = current.size() - start_index;
    
    //double temperature = 1000;
    
    for(int i = 0; i < iterations; ++i) {
        
        do {
            //swap0 = get_random_int(current.size());
            //swap1 = get_random_int(current.size());
            swap0 = get_random_int(total) + start_index;
            swap1 = get_random_int(total) + start_index;
        } while(swap0 == swap1);
        
        // swap random peels
        tmp = current[swap0];
        current[swap0] = current[swap1];
        current[swap1] = tmp;
        
        new_cost = calculate_cost(current);
        
        
        if((new_cost != -1) and (new_cost < cost)) {
            printf("iteration %d: %d\n", i, new_cost);
            iter = i;
        }
        
        
        p.increment();
        
        /*
        if((i % 1000) == 0) {
            temperature *= 0.995;
        }
        
        if((new_cost != -1) and (new_cost < best_cost)) {
            best_cost = new_cost;
        }
        */
        
        // if better, store result
        if((new_cost != -1) and (new_cost <= cost)) {
            cost = new_cost;
            //printf("iteration %d: * %d\n", i, new_cost);
            continue;
        }
        
        /*
        if((new_cost != -1) and (exp((cost - new_cost) / temperature) > get_random())) {
            cost = new_cost;
            printf("iteration %d: %d\n", i, new_cost);
            continue;
        }
        */
        
        // not better swap back
        tmp = current[swap0];
        current[swap0] = current[swap1];
        current[swap1] = tmp;
    }
    
    p.finish_msg("final cost = %d", calculate_cost(current));
    
    //printf("%d %d\n", cost, iter);
    
    //exit(0);
    
    //printf("best = %d\n", best_cost);
    
    //printf("%d\n", calculate_cost(current));
    
    
    for(unsigned int i = 0; i < current.size(); ++i) {
        printf("%d\n", current[i]);
    }
    
    rebuild_peel_order(current);
    
    //exit(0);
    
    
    /*
    // section 6.3, peeling sequence generation
    // greedy
    vector<int> current;
    
    for(unsigned i = 0; i < peelorder.size(); ++i) {
        current.push_back(peelorder[i].get_peelnode());
    }
    
    printf("%d\n", calculate_cost(current));
    
    exit(0);
    */
    
    
    for(unsigned int i = 0; i < peelorder.size(); ++i) {
        bruteforce_assignments(peelorder[i]);
    }
}

bool PeelSequenceGenerator::read_from_file(string filename) {
    SimpleParser sp(filename);
    if(not sp.parse()) {
        fprintf(stderr, "bad parse\n");
        return false;
    }
    
    vector<unsigned int>& current = sp.get_values();
    
    printf("read peeling sequence, cost = %d\n", calculate_cost(current));
    
    rebuild_peel_order(current);
    
    for(unsigned int i = 0; i < peelorder.size(); ++i) {
        bruteforce_assignments(peelorder[i]);
    }
    
    return true;
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
    int total, offset, index;
    
    total = pow(4.0, ndim + 1);
    
    vector<vector<int> > assigns(total, vector<int>(ped->num_members(), -1));
    vector<vector<int> > matrix_indices(map->num_markers());
    vector<vector<int> > presum_indices(map->num_markers());
    vector<int> lod_indices;
    
    vector<unsigned int> cutset(op.get_cutset());
    cutset.push_back(op.get_peelnode());
    
    for(int ind = 0; ind < total; ++ind) {
        index = ind;
        
        for(int i = ndim; i > -1; --i) {
            offset = 1 << (i * 2);
            assigns[ind][cutset[i]] = index / offset;
            index %= offset;
        }
    }
    
    // XXX presum indices, ndim is always one bigger
    for(int locus = 0; locus < int(map->num_markers()); ++locus) {
        for(int i = 0; i < total; ++i) {
            bool valid = true;
            
            for(int j = 0; j < (ndim + 1); ++j) {
                if(not ge.is_legal(cutset[j], locus, assigns[i][cutset[j]])) {
                    valid = false;
                    break;
                }
            }
                        
            if(valid) {
                presum_indices[locus].push_back(i);
            }
        }
    }
    
    
    cutset.pop_back();
    total = pow(4.0, ndim);
    vector<vector<int> > assigns2(total, vector<int>(ped->num_members(), -1));
    
    for(int ind = 0; ind < total; ++ind) {
        index = ind;
        
        for(int i = ndim; i > -1; --i) {
            offset = 1 << (i * 2);
            assigns2[ind][cutset[i]] = index / offset;
            index %= offset;
        }
    }
    
    
    // XXX matrix_indices
    for(int locus = 0; locus < int(map->num_markers()); ++locus) {
        for(int i = 0; i < total; ++i) {
            bool valid = true;
            
            for(int j = 0; j < ndim; ++j) {
                if(not ge.is_legal(cutset[j], locus, assigns2[i][cutset[j]])) {
                    valid = false;
                    break;
                }
            }
                        
            if(valid) {
                matrix_indices[locus].push_back(i);
            }
        }
    }
        
    
    // XXX lod indices
    enum phased_trait pt;
    
    for(int i = 0; i < total; ++i) {
        bool valid = true;
        
        for(int j = 0; j < ndim; ++j) {
            pt = static_cast<enum phased_trait>(assigns[i][cutset[j]]);
            
            if((ped->get_by_index(cutset[j]))->get_disease_prob(pt) == 0.0) {
                valid = false;
                break;
            }
        }
        
        if(valid) {
            lod_indices.push_back(i);
        }
    }
    
    
    op.set_index_values(assigns2);
    op.set_matrix_indices(matrix_indices);
    op.set_presum_indices(presum_indices);
    op.set_lod_indices(lod_indices);
}

int PeelSequenceGenerator::calculate_cost(vector<unsigned int>& seq) {
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

void PeelSequenceGenerator::rebuild_peel_order(vector<unsigned int>& seq) {
    
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

        if(verbose) {
            fprintf(stderr, "rebuild: %s\n", p.debug_string().c_str());
        }
        
        state.toggle_peel_operation(p);
        
        peelorder.push_back(p);
    }
}

