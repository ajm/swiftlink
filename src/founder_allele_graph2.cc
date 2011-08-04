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


bool FounderAlleleGraph::populate_graph(DescentGraph& d) {
    for(unsigned i = 0; i < ped->num_members(); ++i) {
		Person* tmp = ped->get_by_index(i);
        
        if(not tmp->istyped())
			continue;
		
		if(not matrix.add(d.get_founderallele(i, locus, MATERNAL), d.get_founderallele(i, locus, PATERNAL), tmp->get_genotype(locus)))
		    return false;
	}
	
	return true;
}

void FounderAlleleGraph::populate_components() {
    queue<unsigned> q;
    GraphComponent gc;
    vector<int> visited(num_alleles, WHITE);
    int remaining = num_alleles;
    
    do {
        // get a starting point
        for(int i = 0; i < num_alleles; ++i) {
            if(visited[i] == WHITE) {
                visited[i] = GREY;
                q.push(i);
                break;
            }
        }
        
        // breadth first
        while(not q.empty()) {
            int tmp = q.front(); q.pop();
            
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
        
        components.add(gc);
        gc.clear();
        
    } while(remaining != 0);
}

bool FounderAlleleGraph::likelihood(double* prob) {
    double log_prob = 0.0;
    double tmp_prob;
    
    for(unsigned i = 0; i < components.size(); ++i) {
        if((tmp_prob = enumerate_component(component[i])) == 0.0) {
            return false;
        }
        
        log_prob += log(tmp_prob);
    }
    
    return true;
}

double FounderAlleleGraph::component_likelihood(vector<unsigned> q) {
    double tmp = 1.0;
    for(unsigned i = 0; i < q.size(); ++i) {
        tmp *= ((q[i] == 1) ? map->get_major(locus) : map->get_minor(locus));
    }
    return tmp;
}

// for loops in the allele graph, hetero is always a contradiction
bool FounderAlleleGraph::correct_alleles(enum unphased_genotype g, int allele1) {
    
    switch(g) {
        case HOMOZ_A:
            return allele1 == 1;
            
        case HOMOZ_B:
            return allele1 == 2;
        
        case HETERO:
            return false;
    }
}

bool FounderAlleleGraph::correct_alleles(enum unphased_genotype g, int allele1, int allele2) {
    
    switch(g) {
        case HOMOZ_A:
            return (allele1 == 1) and (allele2 == 1);
            
        case HOMOZ_B:
            return (allele1 == 2) and (allele2 == 2);
        
        case HETERO:
            return ((allele1 == 1) and (allele2 == 2)) or \
                   ((allele1 == 2) and (allele2 == 1));
    }
}

// assignment can be of any length
// this code is pretty crufty, maybe warrenting a change of 
// data structures or something...
bool FounderAlleleGraph::legal(GraphComponent& gc, vector<unsigned>& assignment) {
    
    int node = gc[assignment.size() - 1];
    int allele = assignment.back();
    
    AdjacencyRecord tmp = matrix[node];
    
    for(unsigned i = 0; i < tmp.size(); ++i) {
        FounderAlleleNode adj = tmp[i];
        
        // if there is a loop
        if(adj.id == node) {
            if(not correct_alleles(adj.label, allele)) {
                return false;
            }
            
            continue;
        }
        
        // find offset of adjacent node in assignment vector
        unsigned j;
        for(j = 0; j < gc.size(); ++j) {
            if(gc[j] == adj.id)
                break;
        }
        
        // error if not found
        if(j == gc.size()) {
            fprintf(stderr, "Error: an adjacent allele in the graph was not found in the same component (%s:%d)", __FILE__, __LINE__);
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

double FounderAlleleGraph::enumerate_component(GraphComponent& c) {
    vector<unsigned> q;
    double prob = 0.0;
	
	if(c.size() == 1) {
	    return 1.0;
	}

// this generates all permutations, we don't want that...	
/*	
    while(true) {
	    while(q.size() != c.size()) {
	        q.push_back(1);
	    }
	    
	    if(legal(q))
    	    prob += component_likelihood(q);
        
        while((q.size() != 0) and (q.back() == 2)) {
            q.pop();
        }
        
        if(q.size() == 0)
            break;
            
        q.back() = 2;
	}
*/
	
	while(true) {
	    // this while loop does two things, fills the vector with 1's up to the correct length
	    // in the event that a 1 produces an illegal vector, it tries 2 at the same position
	    // if both 1 and 2 are illegal for this allele then there are no valid assignments for 
	    // this component and hence the descent graph is illegal
	    while(q.size() != c.size()) {
	        q.push_back(1);
	        
	        if(not legal(c, q)) {
	            q.back() = 2;
	            if(not legal(c, q)) {
	                goto no_more_assignments;
	            }
	        }
	    }
	    
	    prob += component_likelihood(q);
        
        while(true) {
            // remove all the trailing 2's
            while((q.size() != 0) and (q.back() == 2)) {
                q.pop();
            }
            
            // check to see if we are done
            if(q.size() == 0)
                goto no_more_assignments;
            
            // increment the last allele
            q.back() = 2;
            
            // if this is legal, the we will loop back to the beginning (and fill with 1's)
            // otherwise the current while loop will remove the trailing 2 etc, etc...
            if(legal(c, q))
                break;
        }
	}
    
no_more_assignments:
    
    return prob;
}

