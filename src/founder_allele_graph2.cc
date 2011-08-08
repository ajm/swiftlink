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
    if(not populate_graph(d))
        return false;
        
    populate_components();
    
    return true;
}

bool FounderAlleleGraph2::populate_graph(DescentGraph& d) {
    for(unsigned i = 0; i < ped->num_members(); ++i) {
		Person* tmp = ped->get_by_index(i);
        
        //printf("Adding %d...\n", tmp->get_internalid());
        
        if(not tmp->istyped())
			continue;
		
		if(not matrix.add(d.get_founderallele(i, locus, MATERNAL), d.get_founderallele(i, locus, PATERNAL), tmp->get_genotype(locus)))
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
    
    do {
        // get a starting point
        for(unsigned i = 0; i < num_alleles; ++i) {
            if(visited[i] == WHITE) {
                visited[i] = GREY;
                q.push(i);
                break;
            }
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
            --remaining;
            
            gc.add(tmp);
        }
        
        components.push_back(gc);
        gc.clear();
        
    } while(remaining != 0);
}

bool FounderAlleleGraph2::likelihood(double* prob) {
    double log_prob = 0.0;
    double tmp_prob;
    
    for(unsigned i = 0; i < components.size(); ++i) {
        /*
        for (unsigned j = 0; j < components[i].size(); ++j) {
            fprintf(stderr, "%d ", components[i][j]);
        }
        fprintf(stderr, "\n");
        */
        
        if((tmp_prob = enumerate_component(components[i])) == 0.0) {
            return false;
        }
        
        log_prob += log(tmp_prob);
    }
    
    *prob = log_prob;
    
    return true;
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
        
    AdjacencyRecord tmp = matrix[node];
    
    for(unsigned i = 0; i < tmp.size(); ++i) {
        FounderAlleleNode adj = tmp[i];
                
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

double FounderAlleleGraph2::enumerate_component(GraphComponent& c) {
    vector<unsigned> q;
    double prob = 0.0;
	
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

