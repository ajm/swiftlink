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

using namespace std;


/*
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
*/

vector<PeelOperation>& PeelSequenceGenerator::get_peel_order() {
    return peelorder;
}

string PeelSequenceGenerator::debug_string() {
    stringstream ss;
    
    for(unsigned i = 0; i < peelorder.size(); ++i) {
        ss << i << "\t" << peelorder[i].translated_debug_string(ped) << "\n";
        //ss << peelorder[i].code_output() << "\n";
    }
    
    return ss.str();
}

void PeelSequenceGenerator::find_prev_functions(PeelOperation& op) {
    vector<unsigned> tmp(op.get_cutset());
    tmp.push_back(op.get_peelnode());
    
    while(1) {
        int pf = find_function_containing(tmp);
        
        if(pf == -1)
            break;
            
        op.add_prevfunction(pf);
    }
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
    
    // presum indices, ndim is always one bigger
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
    
    
    // matrix_indices
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
        
    
    // lod indices
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

void PeelSequenceGenerator::set_type(PeelOperation& p) {
    enum peeloperation t = NULL_PEEL;
    
    Person* q = ped->get_by_index(p.get_peelnode());
    
    if(state.is_final_node(p.get_peelnode())) {
        t = LAST_PEEL;
    }
    // not a founder and neither parent has been peeled
    // then; must be child-peel
    else if((not q->isfounder()) and (not (state.is_peeled(q->get_maternalid()) or state.is_peeled(q->get_paternalid())))) {
        t = CHILD_PEEL;
    }
    // if offspring have not been peeled, then their transmission prob can only 
    // be considered if this is a parent peel
    // with the exception of when their other parent has been peeled
    else if(not q->isleaf() and not q->partners_peeled(state) and not q->offspring_peeled(state)) {
        t = PARENT_PEEL;
    }
    else {
        t = PARTNER_PEEL;
    }
    
    /*
    if(t == NULL_PEEL) {
        fprintf(stderr, 
                "Error: invalid peeling sequence (%s did not have a type)\n", 
                q->get_id().c_str());
        abort();
    }
    */
    
    p.set_type(t);
}

// previously peelorder was just in the order of the pedigree, but after finalise_peel_order
// it is in the order of the peeling sequence itself
void PeelSequenceGenerator::finalise_peel_order(vector<unsigned int>& seq) {
    vector<PeelOperation> tmp(peelorder);
    
    peelorder.clear();
    state.reset();

    for(unsigned i = 0; i < seq.size(); ++i) {
        PeelOperation& p = tmp[seq[i]];
        
        set_type(p);
        find_prev_functions(p);
        bruteforce_assignments(p);
        
        eliminate_node(tmp, seq[i]);
        
        peelorder.push_back(p);
        state.toggle_peel_operation(p);
    }
}

void PeelSequenceGenerator::build_simple_graph() {
    
    peelorder.clear();
    
    for(unsigned int i = 0; i < ped->num_members(); ++i) {
        PeelOperation po(ped, i);
        
        Person* p = ped->get_by_index(i);
        
        if(not p->isfounder()) {
            po.add_cutnode(p->get_maternalid());
            po.add_cutnode(p->get_paternalid());
        }
        
        for(unsigned int j = 0; j < p->num_children(); ++j) {
            Person* c = p->get_child(j);
            po.add_cutnode(c->get_internalid());
        }
        
        for(unsigned int j = 0; j < p->num_mates(); ++j) {
            Person* c = p->get_mate(j);
            po.add_cutnode(c->get_internalid());
        }
        
        peelorder.push_back(po);
    }
}

void PeelSequenceGenerator::build_peel_sequence(unsigned int iterations) {
    vector<unsigned int> current;
    
    // attempt to find a greedy solution, this will only return true if the
    // pedigree is outbred and it is therefore optimal
    if(greedy_search(current)) {
        if(is_legit(current)) {
            finalise_peel_order(current);
            return;
        }
    }
    
    // otherwise attempt to optimise a better solution using random
    // downhill searching
    while(1) {
        current.clear();
        current.reserve(ped->num_members());
    
        for(unsigned i = 0; i < ped->num_members(); ++i) {
            current.push_back(i);
        }
        
        random_shuffle(current.begin(), current.end());
        
        random_downhill_search(current, iterations);
        
        if(is_legit(current)) {
            break;
        }
    }
    
    finalise_peel_order(current);
    
    //printf("%s", debug_string().c_str());
}

bool PeelSequenceGenerator::greedy_search(vector<unsigned int>& current) {
    vector<PeelOperation> tmp(peelorder);
    
    while(1) {
        if(current.size() == ped->num_members()) {
            break;
        }
        
        sort(tmp.begin(), tmp.end());
        
        unsigned int cost = tmp[0].get_cutset_size();
        
        if(cost > 2) {
            //printf("[greedy] cancelled\n");
            fprintf(stderr, "Greedy algorithm to find peeling sequence failed, switching to randomised method...\n");
            return false;
        }
        
        vector<PeelOperation>::iterator it = tmp.begin();
        while(it != tmp.end()) {
            if(it->get_cutset_size() != cost) {
                break;
            }
            ++it;
        }
        
        random_shuffle(tmp.begin(), it);
        
        int node = tmp[0].get_peelnode();
        //printf("[greedy] %d (%s)\n", node, ped->get_by_index(node)->get_id().c_str());
        
        //eliminate_node(tmp, node);
        vector<unsigned int>& cutset = tmp[0].get_cutset();
        for(unsigned int i = 1; i < tmp.size(); ++i) {
            if(tmp[i].in_cutset(node)) {
                tmp[i].remove_cutnode(node);
                
                for(unsigned int j = 0; j < cutset.size(); ++j) {
                    if(cutset[j] != tmp[i].get_peelnode()) {
                        tmp[i].add_cutnode(cutset[j]);
                    }
                }
            }
        }
        
        current.push_back(node);
        tmp.erase(tmp.begin());
    }
    
    return true;
}

void PeelSequenceGenerator::random_downhill_search(vector<unsigned int>& current, unsigned int iterations) {
    
    int swap0, swap1, tmp, iter;
    int new_cost, cost;
    
    cost = get_cost(current);
    iter = -1;
    swap0 = swap1 = 0;
    
    
    Progress p("Peel Sequence:", iterations);
    
    for(unsigned int i = 0; i < iterations; ++i) {
        
        do {
            swap0 = get_random_int(current.size());
            swap1 = get_random_int(current.size());
        
        } while(swap0 == swap1);
        
        // swap random peels
        tmp = current[swap0];
        current[swap0] = current[swap1];
        current[swap1] = tmp;
        
        new_cost = get_cost(current);
        //printf("iteration %d: %d\n", i, new_cost);
        
        if(new_cost < cost) {
            //printf("iteration %d: %d\n", i, new_cost);
            iter = i;
        }
        
        p.increment();
        
        // if better, store result
        if(new_cost <= cost) {
            cost = new_cost;
            continue;
        }
        
        // not better, swap back
        tmp = current[swap0];
        current[swap0] = current[swap1];
        current[swap1] = tmp;
    }
    
    p.finish_msg("cost = %d\n", get_proper_cost(current));
}

void PeelSequenceGenerator::eliminate_node(vector<PeelOperation>& tmp, unsigned int node) {
    vector<unsigned int>& cutset = tmp[node].get_cutset();
    
    for(unsigned int i = 0; i < cutset.size(); ++i) {
        tmp[cutset[i]].remove_cutnode(node);
        
        for(unsigned int j = 0; j < cutset.size(); ++j) {
            if(i != j) {
                tmp[cutset[i]].add_cutnode(cutset[j]);
                tmp[cutset[j]].add_cutnode(cutset[i]);
            }
        }
    }
}

// on larger problems starting from a random sequence causes a floating point
// exception during pow, so just use the cutset size as a proxy for the number
// of operations
unsigned int PeelSequenceGenerator::get_cost(vector<unsigned int>& peel) {
    vector<PeelOperation> tmp(peelorder);
    unsigned int cost = 0;
    
    for(unsigned int i = 0; i < peel.size(); ++i) {
        cost += tmp[peel[i]].get_cutset_size();
        eliminate_node(tmp, peel[i]);
    }
    
    return cost;
}

unsigned int PeelSequenceGenerator::get_proper_cost(vector<unsigned int>& peel) {
    vector<PeelOperation> tmp(peelorder);
    unsigned int cost = 0;
    
    for(unsigned int i = 0; i < peel.size(); ++i) {
        cost += tmp[peel[i]].get_cost();
        eliminate_node(tmp, peel[i]);
    }
    
    return cost;
}

bool PeelSequenceGenerator::is_legit(vector<unsigned int>& peel) {
    vector<PeelOperation> tmp(peelorder);
    PeelingState tmpstate(ped);
    
    for(unsigned int i = 0; i < peel.size(); ++i) {
        Person* q = ped->get_by_index(tmp[peel[i]].get_peelnode());
        
        // basically if it could be a CHILD_PEEL (statements 1 & 2)
        // or a PARENT_PEEL (statements 3,4,5)
        // then we have a problem
        // I really doubt this situation would occur in the "optimal"
        // sequence, but let's not take any chances because the downstream
        // code is not designed for it...
        // 
        // please send hate mail to amedlar@gmail.com subject: "you moron"
        if(
            (not q->isfounder()) and 
            (not (tmpstate.is_peeled(q->get_maternalid()) or tmpstate.is_peeled(q->get_paternalid()))) and
            (not q->isleaf()) and 
            (not q->partners_peeled(tmpstate)) and 
            (not q->offspring_peeled(tmpstate))
        ) {
            return false;
        }
        
        eliminate_node(tmp, peel[i]);
        tmpstate.toggle_peel_operation(tmp[peel[i]]);
    }
    
    return true;
}

unsigned int PeelSequenceGenerator::get_peeling_cost() {
    unsigned int cost = 0;
    
    for(unsigned int i = 0; i < peelorder.size(); ++i) {
        cost += peelorder[i].get_cost();
    }
    
    return cost;
}

