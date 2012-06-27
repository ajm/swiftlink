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




void PeelSequenceGenerator::build_peel_order() {
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
}

void PeelSequenceGenerator::find_prev_functions(PeelOperation& op) {
    vector<unsigned> tmp(op.get_cutset());
    tmp.push_back(op.get_peelnode());
    
    op.set_previous_operation(find_function_containing(tmp));
    op.set_previous_operation(find_function_containing(tmp));
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
    }
    
    return cost;
}

void PeelSequenceGenerator::rebuild_peel_order(vector<unsigned int>& seq) {
    vector<SimpleGraph> tmp(graph);
    
    peelorder.clear();
    state.reset();

    for(unsigned i = 0; i < seq.size(); ++i) {
        Person* per = ped->get_by_index(seq[i]);
        PeelOperation p = per->peel_operation(state);
        
        if(p.get_type() == NULL_PEEL) {
            printf("bad rebuild: %s (%s)\n", p.debug_string().c_str(), per->get_id().c_str());
            abort();
        }
        
        vector<unsigned int>& cs = tmp[seq[i]].get_cutset();
        for(unsigned int j = 0; j < cs.size(); ++j) {
            p.add_cutnode(cs[j], per->is_offspring(cs[j]));
        }
        
        //find_previous_functions(p);
        find_prev_functions(p);

        if(verbose) {
            fprintf(stderr, "rebuild: %s\n", p.debug_string().c_str());
        }
        
        printf("%02d: (%s)\n", i, p.translated_debug_string(ped).c_str());
        
        state.toggle_peel_operation(p);
        
        eliminate_node(tmp, seq[i]);
        
        peelorder.push_back(p);
    }
}

// ----

void PeelSequenceGenerator::build_simple_graph() {
    graph.clear();
    
    for(unsigned int i = 0; i < ped->num_members(); ++i) {
        SimpleGraph sg(i);
        
        Person* p = ped->get_by_index(i);
        
        if(not p->isfounder()) {
            sg.add(p->get_maternalid());
            sg.add(p->get_paternalid());
        }
        
        for(unsigned int j = 0; j < p->num_children(); ++j) {
            Person* c = p->get_child(j);
            sg.add(c->get_internalid());
        }
        
        for(unsigned int j = 0; j < p->num_mates(); ++j) {
            Person* c = p->get_mate(j);
            sg.add(c->get_internalid());
        }
        
        //printf("cost of %d = %d\n", i, sg.get_cost());
        
        graph.push_back(sg);
    }
}

void PeelSequenceGenerator::print_graph(vector<SimpleGraph>& g) {
    for(unsigned int i = 0; i < ped->num_members(); ++i) {
        printf("%s\n", g[i].debug_string().c_str());
    }
    printf("\n");
}

void PeelSequenceGenerator::build_peel_sequence() {
    int iterations = 1000;
    vector<unsigned int> current;
    
    /*
    current.reserve(ped->num_members());
    
    for(unsigned i = 0; i < ped->num_members(); ++i) {
        current.push_back(i);
    }
    
    random_shuffle(current.begin(), current.end());
    
    
    int swap0, swap1, tmp, iter;
    int new_cost, cost;
    
    cost = get_cost(current);
    iter = -1;
    swap0 = swap1 = 0;
    
    printf("start cost = %d\n", cost);
    
    //exit(-1);
    
    Progress p("Peel Sequence:", iterations);
    
    for(int i = 0; i < iterations; ++i) {
        
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
            printf("iteration %d: %d\n", i, new_cost);
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
    
    p.finish_msg("final cost = %d (old method = %d)", get_proper_cost(current), calculate_cost(current));
    */ 
    
    /*
    // fixed
    current.clear();
    current.push_back(0);
    current.push_back(32);
    current.push_back(11);
    current.push_back(25);
    current.push_back(14);
    current.push_back(5);
    current.push_back(8);
    current.push_back(15);
    current.push_back(7);
    current.push_back(12);
    current.push_back(23);
    current.push_back(30);
    current.push_back(34);
    current.push_back(28);
    current.push_back(20);
    current.push_back(1);
    current.push_back(13);
    current.push_back(24);
    current.push_back(47);
    current.push_back(43);
    current.push_back(48);
    current.push_back(10);
    current.push_back(39);
    current.push_back(2);
    current.push_back(36);
    current.push_back(46);
    current.push_back(35);
    current.push_back(45);
    current.push_back(44);
    current.push_back(21);
    current.push_back(33);
    current.push_back(3);
    current.push_back(27);
    current.push_back(17);
    current.push_back(4);
    current.push_back(19);
    current.push_back(16);
    current.push_back(6);
    current.push_back(18);
    current.push_back(49);
    current.push_back(50);
    current.push_back(9);
    current.push_back(38);
    current.push_back(37);
    current.push_back(26);
    current.push_back(42);
    current.push_back(31);
    current.push_back(29);
    current.push_back(40);
    current.push_back(41);
    current.push_back(22);
    */
    
    
    // wrong p(T)
    current.clear();
    current.push_back(9);
    current.push_back(7);
    current.push_back(30);
    current.push_back(0);
    current.push_back(2);
    current.push_back(12);
    current.push_back(23);
    current.push_back(6);
    current.push_back(13);
    current.push_back(47);
    current.push_back(16);
    current.push_back(36);
    current.push_back(48);
    current.push_back(46);
    current.push_back(10);
    current.push_back(5);
    current.push_back(15);
    current.push_back(11);
    current.push_back(8);
    current.push_back(14);
    current.push_back(35);
    current.push_back(39);
    current.push_back(44);
    current.push_back(24);
    current.push_back(25);
    current.push_back(21);
    current.push_back(4);
    current.push_back(17);
    current.push_back(19);
    current.push_back(1);
    current.push_back(20);
    current.push_back(49);
    current.push_back(3);
    current.push_back(45);
    current.push_back(27);
    current.push_back(43);
    current.push_back(32);
    current.push_back(28);
    current.push_back(33);
    current.push_back(34);
    current.push_back(18);
    current.push_back(38);
    current.push_back(50);
    current.push_back(22);
    current.push_back(37);
    current.push_back(26);
    current.push_back(29);
    current.push_back(31);
    current.push_back(41);
    current.push_back(42);
    current.push_back(40);
    
    
    /*
    // wrong cost?
    current.clear();
    current.push_back(0);
    current.push_back(24);
    current.push_back(2);
    current.push_back(30);
    current.push_back(7);
    current.push_back(3);
    current.push_back(1);
    current.push_back(33);
    current.push_back(12);
    current.push_back(23);
    current.push_back(8);
    current.push_back(25);
    current.push_back(6);
    current.push_back(4);
    current.push_back(39);
    current.push_back(15);
    current.push_back(48);
    current.push_back(11);
    current.push_back(13);
    current.push_back(14);
    current.push_back(21);
    current.push_back(18);
    current.push_back(44);
    current.push_back(19);
    current.push_back(26);
    current.push_back(20);
    current.push_back(47);
    current.push_back(29);
    current.push_back(43);
    current.push_back(46);
    current.push_back(42);
    current.push_back(45);
    current.push_back(49);
    current.push_back(36);
    current.push_back(27);
    current.push_back(9);
    current.push_back(17);
    current.push_back(10);
    current.push_back(5);
    current.push_back(16);
    current.push_back(22);
    current.push_back(32);
    current.push_back(50);
    current.push_back(37);
    current.push_back(35);
    current.push_back(34);
    current.push_back(28);
    current.push_back(40);
    current.push_back(38);
    current.push_back(31);
    current.push_back(41);
    */
    
    rebuild_peel_order(current);
    
    printf("final cost = %d (old method = %d)\n", get_proper_cost(current), calculate_cost(current));
    
    for(unsigned int i = 0; i < peelorder.size(); ++i) {
        bruteforce_assignments(peelorder[i]);
    }
    
    printf("final cost = %d (old method = %d)\n", get_proper_cost(current), calculate_cost(current));
}

void PeelSequenceGenerator::eliminate_node(vector<SimpleGraph>& tmp, unsigned int node) {
    vector<unsigned int>& cutset = tmp[node].get_cutset();
    
    for(unsigned int i = 0; i < cutset.size(); ++i) {
        tmp[cutset[i]].remove_node(node);
        for(unsigned int j = 0; j < cutset.size(); ++j) {
            if(i != j) {
                tmp[cutset[i]].add(cutset[j]);
                tmp[cutset[j]].add(cutset[i]);
            }
        }
    }
}

// on larger problems starting from a random sequence causes a floating point
// exception during pow
unsigned int PeelSequenceGenerator::get_cost(vector<unsigned int>& peel) {
    vector<SimpleGraph> tmp(graph);
    unsigned int cost = 0;
    
    for(unsigned int i = 0; i < peel.size(); ++i) {
        cost += tmp[peel[i]].get_cutset_size();
        eliminate_node(tmp, peel[i]);
    }
    
    return cost;
}

unsigned int PeelSequenceGenerator::get_proper_cost(vector<unsigned int>& peel) {
    vector<SimpleGraph> tmp(graph);
    unsigned int cost = 0;
    
    for(unsigned int i = 0; i < peel.size(); ++i) {
        printf("cost: %s\n", tmp[peel[i]].debug_string().c_str());
        cost += tmp[peel[i]].get_cost();
        eliminate_node(tmp, peel[i]);
    }
    
    return cost;
}

