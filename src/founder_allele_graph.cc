using namespace std;

#include <cstdio>
#include <cmath>
#include <queue>
#include <vector>
#include <numeric>
#include <algorithm>

#include "lkg.h"
#include "person.h"
#include "pedigree.h"
#include "founder_allele_graph.h"
#include "descent_graph.h"
#include "genetic_map.h"
#include "genotype.h"


FounderAlleleGraph::FounderAlleleGraph(GeneticMap* g, Pedigree* p) : map(g), ped(p) {
    num_founder_alleles = ped->num_founders() * 2;
    num_neighbours = new int[num_founder_alleles];
    adj_matrix = new adj_node*[num_founder_alleles];
    
    for(int i = 0; i < num_founder_alleles; ++i) {
        // +1 to allow for a node to be it's own neighbour like a some
        // cosanguinous pedigrees
        adj_matrix[i] = new adj_node[num_founder_alleles + 1];
    }
    
    best_descentstate = new int[num_founder_alleles];
}

FounderAlleleGraph::~FounderAlleleGraph() {
    for(int i = 0; i < num_founder_alleles; ++i) {
        delete[] adj_matrix[i];
    }
    delete[] adj_matrix;
    delete[] num_neighbours;
    delete[] best_descentstate;
}

void FounderAlleleGraph::reset() {
    for(int i = 0; i < num_founder_alleles; ++i) {
        num_neighbours[i] = 0;
    }
}

int FounderAlleleGraph::_get_num_neighbours(int node) const {
	return num_neighbours[node];
}

bool FounderAlleleGraph::_add(int mat_fa, int pat_fa, enum unphased_genotype g) {
    adj_node* tmp;
    
    if(g == UNTYPED)
        return true;
    
    // check to see if it exists
    // if it does, then the edge needs to have the same label
    for(int i = 0; i < num_neighbours[mat_fa]; ++i) {
        tmp = &(adj_matrix[mat_fa][i]);
        
        if(tmp->id == pat_fa) {
            return tmp->label == g;
        }
    }
    
    tmp = &(adj_matrix[mat_fa][num_neighbours[mat_fa]++]);
    tmp->id = pat_fa;
    tmp->label = g;
    
    tmp = &(adj_matrix[pat_fa][num_neighbours[pat_fa]++]);
    tmp->id = mat_fa;
    tmp->label = g;
        
    return true;
}

bool FounderAlleleGraph::_check_legality(int node, int node_assignment, vector<int>& assignments) {
    adj_node* adj;
    int adj_assignment;
    
    // for each adjacency in the graph
    for(int i = 0; i < num_neighbours[node]; ++i ) {
        adj = &(adj_matrix[node][i]);
        adj_assignment = assignments[adj->id];
        
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

void FounderAlleleGraph::print() const {
    adj_node* tmp;
    
    printf("\nFOUNDER ALLELE GRAPH:\n");
    
	for(int i = 0; i < num_founder_alleles; ++i) {
        printf("%d: ", i);
        for(int j = 0; j < num_neighbours[i]; ++j) {
            tmp = &(adj_matrix[i][j]);
            printf("{%d, %d} ", tmp->id, tmp->label);
        }
        printf("\n");
    }
    printf("\n");
}

bool FounderAlleleGraph::populate(DescentGraph& d, unsigned locus) {
    Person* tmp;    
    int mat_fa, pat_fa ;
	enum unphased_genotype geno;
	
	for(unsigned i = 0; i < ped->num_members(); ++i) {
		tmp = ped->get_by_index(i);
        
        if(not tmp->istyped())
			continue;
		
		mat_fa = d.get_founderallele(i, locus, MATERNAL);
		pat_fa = d.get_founderallele(i, locus, PATERNAL);
		
		geno = tmp->get_genotype(locus);

		//printf("m=%d p=%d l=%u\n", mat_fa, pat_fa, locus);
		
		if(not _add(mat_fa, pat_fa, geno)) {
			return false;
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

	//printf("component_size = %d\n", component_size);
	
    _assign_and_recurse(component, component_size, locus, 
                        0, assignment, &prob, &best);
    
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

double FounderAlleleGraph::descentstate_likelihood(unsigned locus) {
    double prob = 1.0;
    
    for(int i = 0; i < num_founder_alleles; ++i) {
        /*
        prob *= ((best_descentstate[i] == 1) ? \
                    map->get_major(locus) : \
                    map->get_minor(locus));
        */

        switch(best_descentstate[i]) {
            case 1:
                prob *= map->get_major(locus);
                break;
            case 2:
                prob *= map->get_minor(locus);
                break;
            default:
                break;
        }
    }
    
    return log(prob);
}

int FounderAlleleGraph::get_founderallele_assignment(int fa) {
    return best_descentstate[fa];
}

void FounderAlleleGraph::print_ds() {
    for(int i = 0; i < num_founder_alleles; ++i)
        printf("%d ", best_descentstate[i]);
    printf("\n");
}

