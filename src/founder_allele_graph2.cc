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


bool FounderAlleleGraph::populate(DescentGraph& d, unsigned locus) {    
    for(unsigned i = 0; i < ped->num_members(); ++i) {
		Person* tmp = ped->get_by_index(i);
        
        if(not tmp->istyped())
			continue;
		
		if(matrix.add(d.get_founderallele(i, locus, MATERNAL), d.get_founderallele(i, locus, PATERNAL), tmp->get_genotype(locus)))
		    return false;
	}
	
	return true;
}








// this is still old stuff below...




bool FounderAlleleGraph::_check_legality(int node, int node_assignment, vector<int>& assignments) {
    adj_node* adj;
    int adj_assignment;
    
    // for each adjacency in the graph
    for(int i = 0; i < num_neighbours[node]; ++i ) {
        adj = &(adj_matrix[node][i]);
        adj_assignment = assignments[adj->id];
        
        // test for loops in the founder allele graph
        // this will not have been assigned yet
        if(adj->id == node) {
            if(!( \
                ((adj->label == HOMOZ_A) && (node_assignment == 1)) || \
                ((adj->label == HOMOZ_B) && (node_assignment == 2))) \
              ) {
             
                return false;
            }           
        }
        
        if(adj_assignment == -1) { // not yet assigned
            continue;
        }
        
        switch(adj->label) {
            case HOMOZ_A:
                if(adj_assignment != 1 || node_assignment != 1) {
                    return false;
                }
                break;
                
            case HOMOZ_B:
                if(adj_assignment != 2 || node_assignment != 2) {
                    return false;
                }
                break;
                
            case HETERO:
                if((adj_assignment == 1 && node_assignment == 1) || \
                   (adj_assignment == 2 && node_assignment == 2)) {
                    return false;
                }
                break;
                
            default:
                break;
        }
    }
    
    return true;
}




void FounderAlleleGraph::_assign_and_recurse(int *component, int component_size, unsigned locus, 
                                             int current_index, vector<int>& assignment, double *prob, double *best) {
    int allele;
    double tmp = 1.0;
	
    if( current_index == component_size ) {
        // calculate prior probability, include in 'prob'
        for(int i = 0; i < component_size; ++i ) {
            allele = assignment[component[i]];
            tmp *= ((allele == 1) ? map->get_major(locus) : map->get_minor(locus));
        }
        
        // keep track of the best assignment for this component
        if(tmp > *best) {
            *best = tmp;
            for(int i = 0; i < component_size; ++i) {
                best_descentstate[component[i]] = assignment[component[i]];
            }
        }
        
        *prob += tmp;
        return;
    }
	
    for(int i = 1; i <= 2; ++i) { // snps only
        // check that assignment is legal, given other assignments
        if(not _check_legality(component[current_index], i, assignment)) {
            continue;
        }
		        
		assignment[component[current_index]] = i;
        _assign_and_recurse(component, 
							component_size, 
							locus, 
                            current_index + 1, 
							assignment, 
							prob, 
							best);
        assignment[component[current_index]] = -1;
    }
}

double FounderAlleleGraph::_enumerate_component(int *component, int component_size, unsigned locus) {
    //int assignment[num_founder_alleles];
    vector<int> assignment(num_founder_alleles, -1);
    double prob = 0.0;
    double best = 0.0;
    
//    for(int i = 0; i < num_founder_alleles; ++i ) {
//        assignment[i] = -1;
//    }

//	printf("component_size = %d\n", component_size);
//    for(unsigned i = 0; i < component_size; ++i) {
//        printf("%d %d\n", i, component[i]);
//    }
	
	if (component_size == 1)
	    return 1.0;
	
    _assign_and_recurse(component, component_size, locus, 
                        0, assignment, &prob, &best);
    
    /*
    for(unsigned i = 0; i < num_founder_alleles; ++i) {
        printf("%d %d\n", i, assignment[i]);
    }
    printf("\n");
    for(unsigned i = 0; i < component_size; ++i) {
        printf("%d %d\n", i, component[i]);
    }
    */
    
    return prob;
}

bool FounderAlleleGraph::likelihood(double* prob, unsigned locus) {
    queue<unsigned> q;
    int* component = new int[num_founder_alleles];
    int cindex;
    int* visited = new int[num_founder_alleles];
    int tmp;
    int total = num_founder_alleles;
    adj_node* tmp2;
    double tmp_prob;
    double log_prob = 0.0;
    
	// initialise
	for(int i = 0; i < num_founder_alleles; ++i) {
		visited[i] = WHITE;
		best_descentstate[i] = 0;
	}
	
	do {
		cindex = 0;
		
		// find component start point
		//		- if i were not initialised here or used again, this could continue
		//		- from where it left off
		for(int i = 0; i < num_founder_alleles; ++i) {
			if(visited[i] == WHITE) {
				visited[i] = GREY;
                q.push(i);
				break;
			}
		}
		
		// bfs to find component at a time...
		while(not q.empty()) {
            tmp = q.front();
            q.pop();
						
			// look at adjacent nodes for each founder allele node 
			//adjacent to tmp put in Q
			for(int i = 0; i < num_neighbours[tmp]; ++i ) {
				tmp2 = &(adj_matrix[tmp][i]);
				
				if(visited[tmp2->id] == WHITE) {
					visited[tmp2->id] = GREY;
					q.push(tmp2->id);
				}
			}
			
			// colour black when all adjacent nodes are in queue
			// add to component
			visited[tmp] = BLACK;
			total--;
			
			component[cindex++] = tmp;
		}
		
		// if no valid assignment of founder alleles exists, then one 
		// element in a product is zero and hence so is the overall likelihood
		if((tmp_prob = _enumerate_component(component, cindex, locus)) == 0.0) {
		    delete[] visited;
		    delete[] component;
		
			return false;
		}
		
		// product of all legal assignments
		log_prob += log(tmp_prob);
		
	} while(total != 0);
	
	*prob = log_prob;

    delete[] visited;
    delete[] component;
	
	return true;
}

