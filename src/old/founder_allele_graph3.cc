using namespace std;

#include <cstdio>
#include <cmath>
#include <queue>
#include <vector>
#include <numeric>
#include <algorithm>

#include "types.h"
#include "person.h"
#include "pedigree.h"
#include "founder_allele_graph3.h"
#include "descent_graph.h"
#include "genetic_map.h"


// for loops in the allele graph, hetero is always a contradiction
bool FounderAlleleGraph3::correct_alleles_loop(enum unphased_genotype g, int allele1) {
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

bool FounderAlleleGraph3::correct_alleles(enum unphased_genotype g, int allele1) {
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

bool FounderAlleleGraph3::correct_alleles(enum unphased_genotype g, int allele1, int allele2) {
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
bool FounderAlleleGraph3::legal(vector<int>& component, vector<int>& q) {
    int node = component[q.size() - 1];
    int allele = q.back();
    unsigned j;
    
    for(unsigned i = 0; i < num_neighbours[node]; ++i) {
        AdjacentNode& adj = adj_matrix[node][i];
            
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
        for(j = 0; j < component.size(); ++j) {
            if(component[j] == adj.id) {
                break;
            }
        }
            
        // error if not found
        if(j == component.size()) {
            fprintf(stderr, "error in legal()\n");
            abort();
        }
            
        // if not assigned yet, then ignore
        if(j > (q.size() - 1))
            continue;
            
        // if assigned, then test legality
        if(not correct_alleles(adj.label, allele, q[j])) {
                return false;
        }
    }
    
    return true;
}

double FounderAlleleGraph3::component_likelihood(int locus, vector<int>& q) {
    double minor = map->get_minor(locus);
    double major = map->get_major(locus);
    double tmp = 1.0;
    
    for(unsigned i = 0; i < q.size(); ++i) {
        tmp *= ((q[i] == 1) ? major : minor);
    }
    
    return tmp;
}

double FounderAlleleGraph3::enumerate(int locus, vector<int>& component) {
    vector<int> q;
    bool skip = false;
    double prob = 0.0;

    if(component.size() == 1) {
        return 1.0;
    }
    
    while(1) {
        
        while(q.size() != component.size()) {
            q.push_back(1);
            if(not legal(component, q)) {
                skip = true;
                break;
            }
        }
        
        if(not skip) {
            prob += component_likelihood(locus, q);
        }
        
        while(1) {
            // not empty and last value is 2
            while((q.size() != 0) and (q.back() == 2)) { 
                q.pop_back();
            }
            
            if(q.size() == 0) {
                goto no_more_assignments;
            }
            
            q.back() = 2;
            
            if(legal(component, q)) {
                skip = false;
                break;
            }
        }
    }
    
no_more_assignments:
    
    return prob;
}

double FounderAlleleGraph3::likelihood(int locus) {
    queue<int> q;
    vector<int> component;
    vector<bool> visited(founder_alleles, WHITE);
    
    int total = founder_alleles;
    
    double tmp_prob;
    double prob = 0.0; // perform all in log
    //double prob = 1.0;
    
    do {
        component.clear();
        
        // find start point
        for(unsigned i = 0; i < founder_alleles; ++i) {
            if(visited[i] == WHITE) {
                visited[i] = GREY;
                q.push(i);
                break;
            }
        }
            
        while(not q.empty()) {
            int tmp = q.front();
            q.pop();
            
            for(unsigned i = 0; i < num_neighbours[tmp]; ++i) {
                AdjacentNode& tmp2 = adj_matrix[tmp][i];
                    
                if(visited[tmp2.id] == WHITE) {
                    visited[tmp2.id] = GREY;
                    q.push(tmp2.id);
                }
            }
                
            visited[tmp] = BLACK;
            total--;
                
            component.push_back(tmp);
        }
            
        tmp_prob = enumerate(locus, component);
        
        /*
        if(tmp_prob == 0.0) {
            return LOG_ZERO;
        }
        */
        
        tmp_prob = ((tmp_prob == 0.0) ? LOG_ZERO : log(tmp_prob));
        prob = log_product(tmp_prob, prob);
        
        //prob *= tmp_prob;
        
    } while(total != 0);
        
    return prob;
    //return log(prob);
}

void FounderAlleleGraph3::assign(int allele1, int allele2, enum unphased_genotype g) {
    int nn = num_neighbours[allele1];
    AdjacentNode& tmp = adj_matrix[allele1][nn];
    tmp.id = allele2;
    tmp.label = g;
    num_neighbours[allele1] = nn + 1;
}

bool FounderAlleleGraph3::add(int mat_allele, int pat_allele, enum unphased_genotype g) {
    if(g == UNTYPED) {
        return true;
    }
    
    // check to see if it exists
    // if it does, then the edge needs to have the same label
    int neighbours = num_neighbours[mat_allele];
    
    for(int i = 0; i < neighbours; ++i) {
        AdjacentNode& tmp = adj_matrix[mat_allele][i];
        
        if(tmp.id == pat_allele) {
            return tmp.label == g;
        }
    }
    
    assign(mat_allele, pat_allele, g);
    
    if(mat_allele != pat_allele) {
        assign(pat_allele, mat_allele, g);
    }
    
    return true;
}

bool FounderAlleleGraph3::populate(DescentGraph& dg, int locus) {
    Person* p;
    int pid;
    int parent_allele;
    enum unphased_genotype g;
    int mat;
    int pat;
    
	// find founder allele assignments, this is only related to the current 
	// descent graph and not whether people are typed or not
	for(unsigned i = 0; i < ped->num_members(); ++i) {
	    pid = (*sequence)[i];
        p = ped->get_by_index(pid);
	    
	    if(p->isfounder()) {
	        assignment[pid * 2] = pid * 2;
	        assignment[(pid * 2) + 1] = (pid * 2) + 1;
	    }
	    else {
	        parent_allele = dg.get(pid, locus, MATERNAL);
	        assignment[pid * 2] = assignment[ (p->get_maternalid() * 2) + parent_allele ];
	        
	        parent_allele = dg.get(pid, locus, PATERNAL);
	        assignment[(pid * 2) + 1] = assignment[ (p->get_paternalid() * 2) + parent_allele ];
	    }
	}
    
	// construct the actual graph from the assignments and the genotype
	// information
	for(unsigned i = 0; i < ped->num_members(); ++i) {
	    pid = (*sequence)[i];
        p = ped->get_by_index(pid);
	    
	    if(p->istyped()) {
	        mat = assignment[pid * 2];
	        pat = assignment[(pid * 2) + 1];
	        g = p->get_marker(locus);
	        
	        if(not add(mat, pat, g)) {
                return false;
            }
        }
	}
	
	return true;
}

double FounderAlleleGraph3::evaluate(DescentGraph& dg, unsigned locus) {
    for(unsigned i = 0; i < founder_alleles; ++i) {
        num_neighbours[i] = 0;
    }
    
    return populate(dg, locus) ? likelihood(locus) : LOG_ZERO;
}

string FounderAlleleGraph3::debug_string() { 
    stringstream ss;
        
	for(unsigned int i = 0; i < founder_alleles; ++i) {
        ss << i << ": ";
        for(unsigned int j = 0; j < num_neighbours[i]; ++j) {
            AdjacentNode& tmp = adj_matrix[i][j];
            ss << "{" << tmp.id << ", " << tmp.label << "} ";
        }
        ss << "\n";
    }
    //ss << "\n";
    
    return ss.str();
}

