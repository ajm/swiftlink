using namespace std;

#include <cstdio>
#include <cmath>
#include <queue>
#include <vector>
#include <numeric>
#include <algorithm>

#include "misc.h"
#include "person.h"
#include "pedigree.h"
#include "founder_allele_graph2.h"
#include "descent_graph.h"
#include "genetic_map.h"
#include "genotype.h"


bool FounderAlleleGraph2::populate(DescentGraph& d) {
    if(not populate_graph(d)) {
        //fprintf(stderr, "bad populate_graph()\n");
        return false;
    }
    
    populate_components();
        
    return true;
}

bool FounderAlleleGraph2::populate_graph(DescentGraph& d) {
    for(unsigned i = 0; i < ped->num_members(); ++i) {
		Person* tmp = ped->get_by_index(i);
        
        //printf("Adding %d...\n", tmp->get_internalid());
        
        if(not tmp->istyped())
			continue;
		
		if(not matrix.add(d.get_founderallele(i, locus, MATERNAL), 
		                  d.get_founderallele(i, locus, PATERNAL), 
		                  tmp->get_genotype(locus)))
		    return false;
		
		//printf("%s\n\n\n", matrix.debug_string().c_str());
	}
	
	return true;
}

void FounderAlleleGraph2::populate_components() {
    queue<unsigned> q;
    GraphComponent gc;
    vector<int> visited(num_alleles, WHITE);
    int remaining = num_alleles;
    //int starting_node;
    
    do {
        // get a starting point
        for(unsigned i = 0; i < num_alleles; ++i) {
            if(visited[i] == WHITE) {
                visited[i] = GREY;
                //starting_node = i;
                q.push(i);
                break;
            }
        }
        
        // XXX starting_node might not be initialised
        // I could more all the business-code into the
        // for-loop, remove the enclosing do-while-loop
        // and test for successful completion immediately
        // afterwards as well
        //breadth_first_search_component(starting_node, &gc, visited);


        // breadth first
        while(not q.empty()) {
            int tmp = q.front(); 
            q.pop();
            
            // find adjacent nodes
            for(unsigned i = 0; i < matrix[tmp].size(); ++i) {
                FounderAlleleNode f = matrix[tmp][i];
                if(visited[f.id] == WHITE) {
                    visited[f.id] = GREY;
                    q.push(f.id);
                }
            }
            
            // visited, colour black and add to component
            visited[tmp] = BLACK;
            --remaining;
            
            gc.add(tmp);
        }

        //remaining -= gc.size();
        components.push_back(gc);
        gc.clear();
        
    } while(remaining != 0);
}

double FounderAlleleGraph2::likelihood() {
    double tmp_prob;
    
    log_prob = 0.0;
    
    for(unsigned i = 0; i < components.size(); ++i) {
        
        if((tmp_prob = enumerate_component(components[i])) == 0.0) {
            log_prob = LOG_ILLEGAL;
            //fprintf(stderr, "bad component\n");
            return log_prob;
        }
        
        tmp_prob = log(tmp_prob);
        
        components[i].set_prob(tmp_prob);
        
        log_prob += tmp_prob;
    }
    
    return log_prob;
}

double FounderAlleleGraph2::component_likelihood(vector<unsigned>& q) {
    double tmp = 1.0;
    for(unsigned i = 0; i < q.size(); ++i) {
        tmp *= ((q[i] == 1) ? map->get_major(locus) : map->get_minor(locus));
    }
    return tmp;
}

// for loops in the allele graph, hetero is always a contradiction
bool FounderAlleleGraph2::correct_alleles_loop(enum unphased_genotype g, int allele1) {
    
    switch(g) {
        case HOMOZ_A:
            return allele1 == 1;
        case HOMOZ_B:
            return allele1 == 2;
        case HETERO:
            return false;
        case UNTYPED:
            abort();
    }
    
    return false; // to kill warnings
}

bool FounderAlleleGraph2::correct_alleles(enum unphased_genotype g, int allele1) {
    
    switch(g) {
        case HOMOZ_A:
            return allele1 == 1;
        case HOMOZ_B:
            return allele1 == 2;
        case HETERO:
            return true;
        case UNTYPED:
            abort();
    }
    
    return false; // to kill warnings
}


bool FounderAlleleGraph2::correct_alleles(enum unphased_genotype g, int allele1, int allele2) {
    
    switch(g) {
        case HOMOZ_A:
            return (allele1 == 1) and (allele2 == 1);
        case HOMOZ_B:
            return (allele1 == 2) and (allele2 == 2);
        case HETERO:
            return ((allele1 == 1) and (allele2 == 2)) or \
                   ((allele1 == 2) and (allele2 == 1));
        case UNTYPED:
            abort();
    }
    
    return false; // to kill warnings
}

// assignment can be of any length
// this code is pretty crufty, maybe warrenting a change of 
// data structures or something...
bool FounderAlleleGraph2::legal(GraphComponent& gc, vector<unsigned>& assignment) {
    
    int node = gc[assignment.size() - 1];
    int allele = assignment.back();
        
    AdjacencyRecord& tmp = matrix[node];
    
    for(unsigned i = 0; i < tmp.size(); ++i) {
        FounderAlleleNode& adj = tmp[i];
                
        // if there is a loop
        if(adj.id == node) {
            if(not correct_alleles_loop(adj.label, allele)) {
                return false;
            }
            
            continue;
        }
        
        // test if compatible with label
        if(not correct_alleles(adj.label, allele)) {
            return false;
        }
        
        // find offset of adjacent node in assignment vector
        unsigned j;
        for(j = 0; j < gc.size(); ++j) {
            if(gc[j] == adj.id)
                break;
        }
        
        // error if not found
        if(j == gc.size()) {
            fprintf(stderr, "Error: an adjacent allele in the graph was not "
                            "found in the same component (%s:%d)", 
                            __FILE__, __LINE__);
            abort();
        }
        
        // if not assigned yet, then ignore
        if(j > (assignment.size() - 1))
            continue;
        
        // if assigned, then test legality
        if(not correct_alleles(adj.label, allele, assignment[j])) {
            return false;
        }
    }
    
    return true;
}

double FounderAlleleGraph2::enumerate_component(GraphComponent& c) {
    vector<unsigned> q;
    double prob = 0.0;
	
    // * if the size of the component is 1 then the probability is 1.0,
    // * if the component is empty it is because in the reevaluate()
    //   function used up all the founder alleles 
	//if(c.size() < 2) {
	if(c.size() == 1) {
        return 1.0;
	}

    bool skip = false;
    
    while(true) {
        
	    while(q.size() != c.size()) {
	        q.push_back(1);
            if(not legal(c, q)) {
                skip = true;
                break;
            }
	    }
	    
	    if(not skip) {
    	    prob += component_likelihood(q);
        }
        
        while(true) {
            
            while((q.size() != 0) and (q.back() == 2)) {
                q.pop_back();
            }
        
            if(q.size() == 0) {
                goto no_more_assignments;
            }
        
            q.back() = 2;
            
            if(legal(c, q)) {
                skip = false;
                break;
            }
        }
	}
       
no_more_assignments:
    
    return prob;
}

/*
GraphComponent& FounderAlleleGraph2::find_existing_component(int allele) {
    for(unsigned i = 0; i < components.size(); ++i) {
        if(components[i].contains(allele))
            return components[i];
    }
    
    abort();
}

void FounderAlleleGraph2::breadth_first_search_component(int starting_node, GraphComponent& gc, vector<int>& visited) {
    queue<unsigned> q;
    
    if(visited[starting_node] == WHITE) {
        q.push(starting_node);
        visited[starting_node] = GREY;
    }
    
    // breadth first
    while(not q.empty()) {
        int tmp = q.front(); 
        q.pop();
        
        // find adjacent nodes
        for(unsigned i = 0; i < matrix[tmp].size(); ++i) {
            FounderAlleleNode f = matrix[tmp][i];
            if(visited[f.id] == WHITE) {
                visited[f.id] = GREY;
                q.push(f.id);
            }
        }
        
        // visited, colour black and add to component
        visited[tmp] = BLACK;
        gc.add(tmp);
    }
}

// XXX probably this will be more of a reevaluate-and-write so
// at the moment you will call populate() -> likelihood() -> reevaluate()
// to get the conditional distribution of a meiosis
// (to do better you would need to do something like create 2 diff
// objects and apply the one that was sampled, but let's just bear
// this in mind and do it if the M-sampler is such a bottle-neck)
//
// XXX what if when this is called the founder allele graph is illegal? - just abort()
//
// XXX this leaves the founder allele graph object in a bad state
// and can hence be only called once
double FounderAlleleGraph2::reevaluate(DescentGraph& d, unsigned person_id, unsigned locus, enum parentage parent) {
    
    if(log_prob == LOG_ILLEGAL) {
        fprintf(stderr, "error: reevaluate() assumes that the previous founder allele graph was legal (%s:%d)\n", __FILE__, __LINE__);
        abort();
    }
    
    Person* tmp = ped->get_by_index(person_id);
    
    // if not typed then it will not change the founder allele graph
    if(not tmp->istyped()) {
        fprintf(stderr, "reevaluate: not typed\n");
        d.flip_bit(person_id, locus, parent);
        return log_prob;
    }
    
    // needed for the breadth-first search
    vector<int> visited(num_alleles, WHITE);
    
    // get the two current founder alleles
    int fa1 = d.get_founderallele(person_id, locus, parent);
    int fa2 = d.get_founderallele(person_id, locus, parent == MATERNAL ? PATERNAL : MATERNAL);
    
    // find the two graphs that these founder alleles belong to
    // they might be the same component
    GraphComponent& gc1 = find_existing_component(fa1);
    GraphComponent& gc2 = find_existing_component(fa2);
    
    // change the descent graph
    d.flip_bit(person_id, locus, parent);
    
    // find the founder allele that is changed by the flip of the meiosis indicator
    int fa3 = d.get_founderallele(person_id, locus, parent);
    
    fprintf(stderr, ">> old = %d->%d, new = %d->%d\n", fa1, fa2, fa2, fa3);
    // remove the old edge and add the new one
    fprintf(stderr, "BEFORE REMOVE:\n%s\n", matrix.debug_string().c_str());
    matrix.remove(fa1, fa2);
    fprintf(stderr, "AFTER REMOVE:\n%s\n", matrix.debug_string().c_str());
    if(not matrix.add(fa2, fa3, tmp->get_genotype(locus))) {
        fprintf(stderr, "reevaluate() bad populate\n");
        return LOG_ILLEGAL;
    }
    fprintf(stderr, "AFTER ADD:\n%s\n", matrix.debug_string().c_str());
    
    // find the two new components
    // XXX what if they are the same component?
    GraphComponent new_gc1;
    GraphComponent new_gc2;
    
    // if new_gc1 and new_gc2 are the same component, then new_gc2 will be empty
    // which will give a likelihood of 1.0 from enumerate_component(), so 
    // I can pretend that it is always 2 components
    breadth_first_search_component(fa1, new_gc1, visited);
    breadth_first_search_component(fa2, new_gc2, visited);
    
    // calculate likelihood
    double gc1_new_prob = enumerate_component(new_gc1);
    double gc2_new_prob = enumerate_component(new_gc2);
    
    fprintf(stderr, "gc1 = %e, gc2 = %e\n", gc1_new_prob, gc2_new_prob);
    
    // could be illegal
    if((gc1_new_prob == 0.0) or (gc2_new_prob == 0.0)) {
        fprintf(stderr, "reevaluate() bad component\n");
        return LOG_ILLEGAL;
    }
    
    // legal, log everything, calculate the ratio and return the new likelihood
    gc1_new_prob = log(gc1_new_prob);
    gc2_new_prob = log(gc2_new_prob);
    
    double numerator = gc1_new_prob + gc2_new_prob;
    double denominator = gc1.get_prob();
    if(&gc1 != &gc2) {
        denominator += gc2.get_prob();
    }
    
    return log_prob + (numerator - denominator);
}
*/
